from utils import snps_utils, sample_utils as su, config, parse_midas_data
from collections import defaultdict
import sys
import numpy

host_type = sys.argv[1] # within, between

sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
hmp_samples = su.get_sample_names('hmp')
mother_samples = su.get_sample_names('mother')
infant_samples = su.get_sample_names('infant')

cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'olm', 'hmp']
samples = {cohort: su.get_sample_names(cohort) for cohort in cohorts}

num_bootstraps = 1000

min_coverage = config.min_median_coverage
alpha = 0.05 # Confidence interval range for rate estimates
low_divergence_threshold = 2e-04
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs
allowed_variant_types = set(['1D','2D','3D','4D'])

good_species_list = parse_midas_data.load_pickled_good_species_list()

# ===================================================================
# dN/dS info

syn_differences = {}
syn_opportunities = {}
syn_pseudocounts = {}
non_differences = {}
non_pseudocounts = {}
non_opportunities = {}

data = {}
# ===================================================================

for species_name in good_species_list:
	
	sys.stderr.write("\nProcessing %s...\n" % species_name)
	
	# Grab QP samples for this species
	qp_sample_lists = {}
	for cohort in cohorts:
		qp_sample_lists[cohort] = sorted(su.load_qp_samples(samples[cohort], species_name)['qp'])
	
	combined_qp_samples = sorted(su.flatten([qp_sample_lists[cohort] for cohort in cohorts]))
	combined_sample_idx_map = {combined_qp_samples[i] : i for i in range(len(combined_qp_samples))}
	
	# Using all QP samples to threshold on sample size
	if len(combined_qp_samples) < min_sample_size:
		sys.stderr.write("Not enough haploid samples!\n")
		continue	
	
	# Load singleton matrices 
	sys.stderr.write("Loading pre-computed singleton rates for %s...\n" % species_name)
	singleton_rate_map = snps_utils.load_singleton_rate_map(species_name)
	sys.stderr.write("Calculating matrix...\n")		
	
	snp_samples, syn_singleton_matrix, syn_doubleton_matrix, syn_difference_matrix, syn_opportunity_matrix	= snps_utils.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '4D', allowed_samples=combined_qp_samples)
	snp_samples, non_singleton_matrix, non_doubleton_matrix, non_difference_matrix, non_opportunity_matrix	= snps_utils.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '1D', allowed_samples=combined_qp_samples)
	
	singleton_matrix = syn_singleton_matrix + non_singleton_matrix
	doubleton_matrix = syn_doubleton_matrix + non_doubleton_matrix
	difference_matrix = syn_difference_matrix + non_difference_matrix
	opportunity_matrix = syn_opportunity_matrix + non_opportunity_matrix
	
	substitution_rate_matrix = difference_matrix*1.0/(opportunity_matrix+(opportunity_matrix==0))
	syn_substitution_rate_matrix = syn_difference_matrix*1.0/(syn_opportunity_matrix+(syn_opportunity_matrix==0))
	
	sys.stderr.write("Done!\n")
	
	for i in range(0, syn_difference_matrix.shape[0]):
		for j in range(i+1, syn_difference_matrix.shape[0]):						
			if syn_opportunity_matrix[i,j]>0 and non_opportunity_matrix[i,j]>0:
				sample_i, sample_j = snp_samples[i], snp_samples[j]
				tp_pair = su.sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples)
				
				if host_type == 'within':
					if sample_subject_map[sample_i] != sample_subject_map[sample_j]:
						continue
				elif host_type == 'between':
					if sample_subject_map[sample_i] == sample_subject_map[sample_j]:
						continue
				
				for dnds_dict in [syn_differences, syn_pseudocounts, syn_opportunities, non_differences, non_pseudocounts, non_opportunities]:
					if tp_pair not in dnds_dict:
						dnds_dict[tp_pair] = defaultdict(list)
				
				syn_differences[tp_pair][species_name].append(syn_difference_matrix[i,j]+1)
				syn_pseudocounts[tp_pair][species_name].append(syn_opportunity_matrix[i,j]*1.0/(syn_opportunity_matrix[i,j]+non_opportunity_matrix[i,j]))
				syn_opportunities[tp_pair][species_name].append(syn_opportunity_matrix[i,j])
				
				non_differences[tp_pair][species_name].append(non_difference_matrix[i,j] )
				non_pseudocounts[tp_pair][species_name].append(non_opportunity_matrix[i,j]*1.0/(syn_opportunity_matrix[i,j]+non_opportunity_matrix[i,j]) )
				non_opportunities[tp_pair][species_name].append(non_opportunity_matrix[i,j])

for dnds_dict in [syn_differences, syn_pseudocounts, syn_opportunities, non_differences, non_pseudocounts, non_opportunities]:
		for tp_pair in dnds_dict:
			for species_name in dnds_dict[tp_pair]:
				dnds_dict[tp_pair][species_name] = numpy.array(dnds_dict[tp_pair][species_name])

import pickle
ddir = config.data_directory
pdir = "%s/pickles" % ddir

dnds_dict_names = {'syn_differences': syn_differences, 'syn_pseudocounts': syn_pseudocounts, 'syn_opportunities': syn_opportunities, 'non_differences': non_differences, 'non_pseudocounts': non_pseudocounts, 'non_opportunities': non_opportunities}

for name in dnds_dict_names:
	pickle.dump(dnds_dict_names[name], open("%s/dnds/%s_%s.pkl" % (pdir, host_type, name), 'wb'))
