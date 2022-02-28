import matplotlib	 
matplotlib.use('Agg') 
from utils import sample_utils as su, config, parse_midas_data
import os.path, sys, numpy
import pylab
from utils import diversity_utils, gene_diversity_utils, temporal_changes_utils, substitution_rates_utils, stats_utils, sfs_utils
from collections import defaultdict
from numpy.random import choice, binomial

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, random, choice, multinomial, shuffle
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.patches import Patch

import pickle
from os import path

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--prevalence_cohort", help="Cohort to base prevalence calculation on", default="all")
parser.add_argument("--within_host_type", help="Timepoint pairs to choose for infant-infant", default="nonconsecutive")
parser.add_argument("--within_mi_pair_type", help="Timepoint pairs to choose for mother-infant", default="all")
parser.add_argument("--sweep_type", help="Full or partial sweep", default="full")
parser.add_argument("--modification_threshold", help="Threshold number of SNP changes to consider modification", default=config.modification_difference_threshold)
args = parser.parse_args()

# Cohort to base prevalence calc on: all, infant, hmp
prevalence_cohort = args.prevalence_cohort
if prevalence_cohort not in ['all', 'infant', 'hmp', 'mother', 'premie', 'nonpremie']:
	sys.exit("Invalid prevalence_cohort. Choose from all, infant, hmp, mother, premie, nonpremie")

# Timepoint pairs for infant-infant: consecutive, longest, nonconsecutive, longest_with_data
within_host_type = args.within_host_type
if within_host_type not in ['consecutive', 'longest', 'nonconsecutive', 'longest_with_data']:
	sys.exit("Invalid within_host_type. Choose from consecutive, longest, nonconsecutive")

# Timepoint pairs for mother-infant: first, last, random, all
within_mi_pair_type = args.within_mi_pair_type
if within_mi_pair_type not in ['first', 'last', 'random', 'all']:
	sys.exit("Invalid within_mi_pair_type. Choose from first, last, random, all")

# Full (0.6) or partial (0.3) sweep: full, partial
sweep_type = args.sweep_type
if sweep_type not in ['full', 'partial']:
	sys.exit("Invalid sweep_type. Choose from full, partial")
if sweep_type == 'full':
	lower_threshold, upper_threshold = 0.2, 0.8
elif sweep_type == 'partial':
	lower_threshold, upper_threshold = 0.35, 0.65

# Modification threshold
modification_difference_threshold = int(args.modification_threshold)
if modification_difference_threshold < 1:
	sys.exit("Invalid modification_threshold. Must be at least 1.")

################################################################################

######
#
# Helper function for calculating the prevalence of the sweeping allele in the larger cohort
#
######
def get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map):
		gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
		
		f1 = A1*1.0/D2
		f2 = A2*1.0/D2
		
		is_reversion = (f1>f2)
		location_tuple = (contig, position)
		is_private_snv = (location_tuple in private_snv_map)
		
		# Now calculate frequency-stratified version
		
		if location_tuple in snv_freq_map:
				f = snv_freq_map[location_tuple]
		else:
				sys.stderr.write("SNP not in map. Shouldn't happen!\n")
				f = -0.5
		
		# Let's impose that private snvs have zero freq (specifically, lim 0^-)				 
		if is_private_snv:
				f = -0.5
		
		# Change f so that it represents
		# frequency of allele at second timepoint
		if is_reversion:
				f = 1-f
		
		return f

# Label format function
def format_tp_label(cohort, tp_pair):
	tpa, tpb = tp_pair
	if tpa[0] == 'M' and tpb[0] == 'I' or tpa[0] == 'I' and tpb[0] == 'M':
		tp1, tp2 = (tpa, tpb) if tpa[0] == 'M' else (tpb, tpa)
	else:
		tp1, tp2 = (tpa, tpb) if (tpa < tpb) else (tpb, tpa)
	subj1, subj2 = tp1[0], tp2[0]
	order1, order2 = int(tp1[1]), int(tp2[1])
	if cohort == 'backhed':
		ilabels = ['I_Birth', 'I_4mon','I_12mon']
		mlabels = ['M_Birth']
	elif cohort == 'ferretti':
		ilabels = ['I_1day','I_3days', 'I_1wk', 'I_1mon', 'I_4mon']
		mlabels = ['M_Birth']
	elif cohort == 'yassour':
		ilabels = ['I_Birth', 'I_2wk','I_1mon', 'I_2mon', 'I_3mon']
		mlabels = ['M_Gest', 'M_Birth', 'M_3mon']
	elif cohort == 'hmp':
		return (str(order1) + " > " + str(order2))
	else:
		return 'error: bad cohort'
	label1 = mlabels[order1-1] if subj1 == 'M' else ilabels[order1-1]
	label2 = mlabels[order2-1] if subj2 == 'M' else ilabels[order2-1]
	return (label1.replace('_','-') + " > " + label2.replace('_','-'))

# Convert timepoint pair frozenset to pair type name
def timepoint_pair_type(tp_pair):
	tp_pair_name_map = {'mother-mother': 'MM', 'infant-infant': 'II', 'mother-infant': 'MI', 'adult-adult': 'AA'}
	
	tp1, tp2 = tp_pair
	s1, s2 = tp1[0], tp2[0]
	
	if s1 == s2:
		return (s1 + s2) # MM, II or AA
	elif (s1 + s2) == 'MI' or (s1 + s2) == 'IM':
		return 'MI'
	else: return 'Invalid timepoint pair'

# Parameters
mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']	= False
mpl.rcParams['legend.fontsize']	 = 'small'

min_sample_size = 3
variant_types = ['1D','4D']
min_snp_change_sample_size = 5

# modification_difference_threshold = config.modification_difference_threshold
replacement_difference_threshold = config.replacement_difference_threshold

# For partitioning SNVs according to prevalence
derived_freq_bins = numpy.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = numpy.arange(0,len(derived_freq_bins)-1)
derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
derived_virtual_xticklabels = ['0','.01','.1','.5','.9','.99','1']

# For partitioning genes into different prevalence classes
gene_freq_bins = numpy.array([-1,0.1,0.5,0.9,2])
gene_freq_xticks      = [-4, -3,  -2,   -1,   0,   1,    2,   3, 4]
gene_freq_xticklabels = ['0','0.1','0.5', '0.9','1','0.9','0.5', '0.1','0']
gene_gain_virtual_freqs = numpy.array([3.5,2.5,1.5,0.5])
gene_loss_virtual_freqs = numpy.array([-3.5,-2.5,-1.5,-0.5])

# Cohorts
cohorts = ['backhed', 'ferretti', 'yassour', 'hmp']
mi_cohorts = ['backhed', 'ferretti', 'yassour']

# Timepoint pair types
tp_pair_names = ['MM', 'MI', 'II', 'AA']

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Pickled QP data
all_samples_pickle_fn = '%s/pickles/all_samples.pkl' % config.data_directory
all_samples = pickle.load(open(all_samples_pickle_fn, 'rb'))

qp_samples_pickle_fn = '%s/pickles/qp_samples.pkl' % config.data_directory
qp_samples = pickle.load(open(qp_samples_pickle_fn, 'rb'))

# HMP samples
hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

# All mother-infant samples
all_mi_samples = su.get_sample_names('all','all')
mi_samples = {cohort: su.get_sample_names(cohort, 'all') for cohort in mi_cohorts}
mother_samples = su.get_sample_names('mother','all')
infant_samples = su.get_sample_names('infant','all')

# Species list
good_species_list = parse_midas_data.parse_good_species_list()

# Species SNP/gene change distributions
# Nested dictionaries: timepoint pair followed by species
species_snp_change_distribution = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}
species_snp_nerrs = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}
species_gene_change_distribution = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}
species_gene_nerrs = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}

species_snp_changes = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}

# Pooled SNP/gene change distributions
pooled_snp_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
pooled_gene_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
pooled_snp_length_distribution = {cohort: defaultdict(list) for cohort in cohorts}
pooled_gene_length_distribution = {cohort: defaultdict(list) for cohort in cohorts}

# typical value, median other sample
pooled_between_snp_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
pooled_between_gene_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}

# closest other sample
pooled_min_between_snp_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
pooled_min_between_gene_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}

# Sums up number of SNP changes based on variant type
total_freq_snps = {tp_pair: {cohort: {} for cohort in cohorts} for tp_pair in tp_pair_names}
total_null_freq_snps = {tp_pair: {cohort: {} for cohort in cohorts} for tp_pair in tp_pair_names}
# Sums up number of SNP changes based on timepoint pair type
pooled_snp_change_by_tp_pair = {cohort: {} for cohort in cohorts}

for tp_pair in tp_pair_names:
	for cohort in cohorts:
			# This holds the # of SNVs in each prevalence class of each var type
			total_freq_snps[tp_pair][cohort] = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}
			
			# This holds the expectation of the # of SNVs in each prevalence class of 1D and 4D
			# (i.e., relative fraction of opportunities for different species), conditioned on prevalence
			total_null_freq_snps[tp_pair][cohort] = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}

# Genes which undergo SNP modifications
# Combine all "cohorts"
modified_gene_subjects = defaultdict(set)

# This sums up across the different var_type categories
total_freq_all_snps = {tp_pair: {cohort: numpy.zeros_like(derived_virtual_freqs) for cohort in cohorts} for tp_pair in tp_pair_names}

# This is the null distribution of prevalence (conditioned on total # of SNVs)
total_null_freq_all_snps = {tp_pair: {cohort: numpy.zeros_like(derived_virtual_freqs)*1.0 for cohort in cohorts} for tp_pair in tp_pair_names}

# Gene gains and losses (only losses under de novo mutation assumption)
total_freq_gains = {tp_pair: {cohort: numpy.zeros(len(gene_freq_bins)-1)*1.0 for cohort in cohorts} for tp_pair in tp_pair_names}
total_freq_losses = {tp_pair: {cohort: numpy.zeros_like(total_freq_gains[tp_pair][cohort]) for cohort in cohorts} for tp_pair in tp_pair_names}
total_null_freq_losses = {tp_pair: {cohort: numpy.zeros_like(total_freq_gains[tp_pair][cohort]) for cohort in cohorts} for tp_pair in tp_pair_names}

# SNV prevalences resolved by subject
# so that we can do bootstrapping	
snv_prevalence_count_map = {tp_pair: {cohort: {} for cohort in cohorts} for tp_pair in tp_pair_names}
snv_prevalence_map = {tp_pair: {cohort: {} for cohort in cohorts} for tp_pair in tp_pair_names}
variant_type_prevalence_map = {tp_pair: {cohort: {'1D':[], '4D':[]} for cohort in cohorts} for tp_pair in tp_pair_names}
# cohort_output_strs = {tp_pair: {cohort : [] for cohort in cohorts} for tp_pair in tp_pair_names}

# dN/dS by cohort and timepoint pair type
dNdS_by_tp_pair = {cohort: defaultdict(list) for cohort in cohorts}

# dN/dS by cohort and prevalence (not of samples, but for groups of SNPs)
dNdS_vars = ['N', 'S', 'N_opp', 'S_opp']
dNdS_by_prev = {tp_pair: {cohort: {var: numpy.zeros_like(derived_virtual_freqs) for var in dNdS_vars} for cohort in cohorts} for tp_pair in tp_pair_names}

# Pickles
total_freq_snps_pickle_fn = '%s/pickles/total_freq_snps.pkl' % config.data_directory
total_freq_all_snps_pickle_fn = '%s/pickles/total_freq_all_snps.pkl' % config.data_directory
total_null_freq_snps_pickle_fn = '%s/pickles/total_null_freq_snps.pkl' % config.data_directory
total_null_freq_all_snps_pickle_fn = '%s/pickles/total_null_freq_all_snps.pkl' % config.data_directory
pooled_snp_pickle_fn = '%s/pickles/pooled_snp_change.pkl' % config.data_directory
pooled_gene_pickle_fn = '%s/pickles/pooled_gene_change.pkl' % config.data_directory
dNdS_by_prev_pickle_fn = '%s/pickles/dNdS_by_prev.pkl' % config.data_directory
snp_change_by_tp_pair_fn = '%s/pickles/snp_change_by_tp_pair.pkl' % config.data_directory
total_freq_gains_fn = '%s/pickles/total_freq_gains.pkl' % config.data_directory
total_freq_losses_fn = '%s/pickles/total_freq_losses.pkl' % config.data_directory
total_null_freq_losses_fn = '%s/pickles/total_null_freq_losses.pkl' % config.data_directory

if False:
	total_freq_snps = pickle.load(open(total_freq_snps_pickle_fn, 'rb'))
	total_freq_all_snps = pickle.load(open(total_freq_all_snps_pickle_fn, 'rb'))
	total_null_freq_snps = pickle.load(open(total_null_freq_snps_pickle_fn, 'rb'))
	total_null_freq_all_snps = pickle.load(open(total_null_freq_all_snps_pickle_fn, 'rb'))
	pooled_snp_change_distribution = pickle.load(open(pooled_snp_pickle_fn, 'rb'))
	pooled_gene_change_distribution = pickle.load(open(pooled_gene_pickle_fn, 'rb'))
	dNdS_by_prev = pickle.load(open(dNdS_by_prev_pickle_fn, 'rb'))
	pooled_snp_change_by_tp_pair = pickle.load(open(snp_change_by_tp_pair_fn, 'rb'))
	total_freq_gains = pickle.load(open(total_freq_gains_fn, 'rb'))
	total_freq_losses = pickle.load(open(total_freq_losses_fn, 'rb'))
	total_null_freq_losses = pickle.load(open(total_null_freq_losses_fn, 'rb'))
else:
	for species_name in good_species_list:
		sys.stderr.write("\nProcessing %s...\n" % species_name)	
		qp_sample_sets = {cohort: [all_samples[index] for index in qp_samples[cohort][species_name]] for cohort in cohorts}
		qp_sample_lists = {cohort: sorted(qp_sample_sets[cohort]) for cohort in cohorts}
		combined_qp_samples = sorted(su.flatten([qp_sample_lists[cohort] for cohort in cohorts]))
		combined_sample_idx_map = {combined_qp_samples[i] : i for i in xrange(0,len(combined_qp_samples))}
		
		# Using Backhed to threshold on sample size
		sample_size = len(qp_sample_sets['backhed'])
		if sample_size < min_sample_size:
			continue
		
		import calculate_private_snvs
		private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)
		
		import calculate_snp_prevalences
		snv_freq_map = calculate_snp_prevalences.parse_population_freqs(prevalence_cohort,species_name,polarize_by_consensus=True)
		snv_freq_keys = snv_freq_map.keys()
		snv_freq_values = snv_freq_map.values()
		
		import core_gene_utils
		gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
		gene_freq_values = numpy.array(gene_freq_map.values())
		gene_freq_weights = gene_freq_values*1.0/gene_freq_values.sum()
		
		# Combined Backhed, Ferretti, Yassour calculations
		sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
		substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
		if substitution_rate_map == {}: # Not enough haploid samples
			sys.stderr.write("Not enough haploid samples!\n")
			continue
		sys.stderr.write("Calculating SNV matrix...\n")
		dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=combined_qp_samples)
		
		snp_difference_matrix = snp_mut_difference_matrix+snp_rev_difference_matrix
		snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
		snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading gene matrix...\n")
		gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=combined_qp_samples)
		gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
		gene_opportunity_matrix = gene_loss_opportunity_matrix
		gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading 1D & 4D opportunity matrices...\n")
		
		difference_matrices = {}
		opportunity_matrices = {}
		
		for var_type in variant_types:		
			dummy_samples, difference_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, var_type, allowed_samples=combined_qp_samples)
			
			difference_matrices[var_type] = difference_matrix
			opportunity_matrices[var_type] = opportunity_matrix
		
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading pre-computed temporal changes...\n")
		temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
		sys.stderr.write("Done!\n")
		
		### Now loop over different cohorts
		for cohort in cohorts:		
			desired_samples = qp_sample_lists[cohort]
			
			within_host_type_param = 'nonconsecutive' if (within_host_type == 'longest_with_data') else within_host_type
			
			same_subject_idxs = su.calculate_mi_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type_param, one_per_mi_pair=False)
			
			#apply_sample_index_map_to_indices(sample_idx_map, idxs):
			#new_idxs = (numpy.array([sample_idx_map[i] for i in idxs[0]]), numpy.array([sample_idx_map[i] for i in idxs[1]]))
			
			# Filter out MI and II timepoints for one per host
			if within_host_type == 'longest_with_data' and (within_mi_pair_type == 'first' or within_mi_pair_type == 'last'):
				idxs1, idxs2 = [], []
				
				# Stores "length" of farthest-apart timepoint pair for each host
				longest_host_tp_pair = {}
				# Same but for mother-infant
				best_mi_host_tp_pair = {}
				
				for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
					idx_i = same_subject_idxs[0][sample_pair_idx]
					idx_j = same_subject_idxs[1][sample_pair_idx]
					sample_i = desired_samples[idx_i]
					sample_j = desired_samples[idx_j]
					subject = sample_order_map[sample_i][0][:-2]				
					tp_pair_length = sample_order_map[sample_j][1] - sample_order_map[sample_i][1]
					
					L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
					
					if L<config.min_opportunities:
							continue
					
					num_snp_changes = len(mutations)+len(reversions)					
					
					if num_snp_changes <= modification_difference_threshold:
						if within_host_type == 'longest_with_data' and (sample_i in infant_samples and sample_j in infant_samples):
							if (subject not in longest_host_tp_pair) or (tp_pair_length > longest_host_tp_pair[subject][1]):
								longest_host_tp_pair[subject] = ((idx_i, idx_j), tp_pair_length)
						elif (sample_i in mother_samples and sample_j in infant_samples):
							# Only consider Yassour MI pairs with mother-delivery
							if cohort == 'yassour' and sample_order_map[sample_i][1] != 2:
								continue
							if within_mi_pair_type == 'first':
								if (subject not in best_mi_host_tp_pair) or (tp_pair_length < best_mi_host_tp_pair[subject][1]):
									best_mi_host_tp_pair[subject] = ((idx_i, idx_j), tp_pair_length)
							elif within_mi_pair_type == 'last':
								if (subject not in best_mi_host_tp_pair) or (tp_pair_length > best_mi_host_tp_pair[subject][1]):
									best_mi_host_tp_pair[subject] = ((idx_i, idx_j), tp_pair_length)
					
				for host in longest_host_tp_pair:
					idxs1.append(longest_host_tp_pair[host][0][0])
					idxs2.append(longest_host_tp_pair[host][0][1])
				for host in best_mi_host_tp_pair:
					idxs1.append(best_mi_host_tp_pair[host][0][0])
					idxs2.append(best_mi_host_tp_pair[host][0][1])
				same_subject_idxs = (numpy.array(idxs1,dtype=numpy.int32), numpy.array(idxs2,dtype=numpy.int32))
			
			# Loop over different pairs of within-host samples
			for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):			 
				
				sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
				sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
				
				subject = sample_order_map[sample_i][0][:-2]
				
				if sample_i and sample_j in hmp_samples:
					tp_i, tp_j = ('A' + str(sample_order_map[sample_i][1]), 'A' + str(sample_order_map[sample_j][1]))
				else:
					tp_i = ('M' if sample_i in mother_samples else 'I') + str(sample_order_map[sample_i][1])
					tp_j = ('M' if sample_j in mother_samples else 'I') + str(sample_order_map[sample_j][1])
				
				tp_pair = frozenset((tp_i, tp_j))
				tp_pair_type = timepoint_pair_type(tp_pair)
				
				i = combined_sample_idx_map[sample_i]
				j = combined_sample_idx_map[sample_j]
				
				# print(sample_order_map[sample_i][0] + " " + str(sample_order_map[sample_i][1]) + " > " + sample_order_map[sample_j][0] + " " + str(sample_order_map[sample_j][1]))
				
				good_idxs = su.calculate_samples_in_different_subjects(subject_sample_map, combined_qp_samples, sample_i)
				good_idxs *= ((snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5))
				
				if good_idxs.sum() < 1:
						continue
				
				L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
				
				if L<config.min_opportunities:
						continue
				
				nerr = L*perr
				
				num_mutations = len(mutations)
				num_reversions = len(reversions)
				num_snp_changes = num_mutations+num_reversions
				
				gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
				
				gene_nerr = gene_L*gene_perr
				num_gains = len(gains)
				num_losses = len(losses)
				num_gene_changes = num_gains+num_losses
				
				if (perr<-0.5) or (gene_perr < -0.5):
						continue
				
				if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
						continue # Only take things with low-ish FPR
				
				# Species specific distributions
				species_snp_change_distribution[cohort][tp_pair][species_name].append( num_snp_changes)
				species_snp_nerrs[cohort][tp_pair][species_name].append(nerr)
				
				species_gene_change_distribution[cohort][tp_pair][species_name].append(num_gene_changes)
				species_gene_nerrs[cohort][tp_pair][species_name].append(gene_nerr)
				
				# Pooled distributions
				pooled_snp_change_distribution[cohort][tp_pair].append(num_snp_changes)
				pooled_gene_change_distribution[cohort][tp_pair].append(num_gene_changes)
				
				pooled_snp_length_distribution[cohort][tp_pair].append(L)
				pooled_gene_length_distribution[cohort][tp_pair].append(gene_L)
				
				# Matched between-host samples
				# typical
				pooled_between_snp_change_distribution[cohort][tp_pair].append( choice(snp_difference_matrix[i, good_idxs]) )
				pooled_between_gene_change_distribution[cohort][tp_pair].append( choice(gene_difference_matrix[i, good_idxs]) )
				
				# minimum
				pooled_min_between_snp_change_distribution[cohort][tp_pair].append( snp_difference_matrix[i, good_idxs].min() )
				pooled_min_between_gene_change_distribution[cohort][tp_pair].append( gene_difference_matrix[i, good_idxs].min() )
				
				subject_type_1 = tp_i[0]
				subject_type_2 = tp_j[0]
				# sys.stderr.write("1:%s 2:%s\n" % (subject_type_1, subject_type_2))
				
				nonsyn_diff_count = 0
				syn_diff_count = 0
				
				# Skip mother-gestation, 3mo postpartum
				if cohort == 'yassour' and ('M1' in tp_pair or 'M3' in tp_pair):
					continue
				
				# If deemed a modification, investigate properties of SNVs and genes
				if (num_snp_changes<=modification_difference_threshold):
						
						# Compute non/syn opportunities for this sample pair
						total_non_opportunities = opportunity_matrices['1D'][i,j]
						total_syn_opportunities = opportunity_matrices['4D'][i,j]
						
						freq_list = set()
						
						for snp_change in (mutations+reversions):
								variant_type = snp_change[3]							
								gene = snp_change[0]
								sample_pair = frozenset((sample_i, sample_j))
								# Previously just added subject, but lost too much info
								modified_gene_subjects[gene].add(sample_pair)
								
								f = get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map)
								
								f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()	
								freq_list.add(f_idx)
								
								if tp_pair not in pooled_snp_change_by_tp_pair[cohort]:
										pooled_snp_change_by_tp_pair[cohort][tp_pair] = numpy.zeros_like(derived_virtual_freqs)
								pooled_snp_change_by_tp_pair[cohort][tp_pair][f_idx] += 1
								
								if variant_type in variant_types:
										total_freq_snps[tp_pair_type][cohort][variant_type][f_idx] += 1
										
										# Calculate null version based on # of opportunities
										total_opportunities = 0.0
										for other_variant_type in variant_types:
												total_opportunities = opportunity_matrices[other_variant_type][i,j]
										
										for other_variant_type in variant_types:
												total_null_freq_snps[tp_pair_type][cohort][other_variant_type][f_idx] += opportunity_matrices[other_variant_type][i,j]/total_opportunities
															
								total_freq_all_snps[tp_pair_type][cohort][f_idx] += 1
								
								if (sample_i, sample_j, species_name) not in snv_prevalence_count_map[tp_pair_type][cohort]:
										snv_prevalence_count_map[tp_pair_type][cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_all_snps[tp_pair_type][cohort])
										snv_prevalence_map[tp_pair_type][cohort][(sample_i, sample_j, species_name)] = []
								
								snv_prevalence_count_map[tp_pair_type][cohort][(sample_i, sample_j, species_name)][f_idx] += 1
								
								snv_prevalence_map[tp_pair_type][cohort][(sample_i, sample_j, species_name)].append(f) 
								
								if variant_type not in variant_type_prevalence_map[tp_pair_type][cohort]:
										variant_type_prevalence_map[tp_pair_type][cohort][variant_type] = []
								
								variant_type_prevalence_map[tp_pair_type][cohort][variant_type].append(f)
								
								# Further investigate dN/dS of SNPs involved in putative modification							
								if variant_type == '1D':
									nonsyn_diff_count += 1
									dNdS_by_prev[tp_pair_type][cohort]['N'][f_idx] += 1
								elif variant_type == '4D':
									syn_diff_count += 1
									dNdS_by_prev[tp_pair_type][cohort]['S'][f_idx] += 1
								
								# Now draw a null prevalence from the genome
								L = snp_opportunity_matrix[i,j]
								L_snv = len(snv_freq_map) # A slight overestimate
								snv_fraction = L_snv*1.0/L
								num_bootstraps = 10
								for bootstrap_idx in xrange(0,num_bootstraps):
										
										if random()<snv_fraction:
												# A polymorphic site
												
												random_snv_idx = randint(0,len(snv_freq_keys))
												random_snv_location = snv_freq_keys[random_snv_idx]
												f = snv_freq_values[random_snv_idx]
												
												rev_f = 1-f
												
												if random_snv_location in private_snv_map:
														# A private SNV. Use private bins
														# use private					
														f_idx = 0
														rev_f_idx = -1
												else:													
														f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()
														rev_f_idx = ((rev_f>derived_freq_bins[:-1])*(rev_f<=derived_freq_bins[1:])).argmax()
												
												# Now add in probability weight
												total_null_freq_all_snps[tp_pair_type][cohort][f_idx] += (1-f)*1.0/num_bootstraps
												total_null_freq_all_snps[tp_pair_type][cohort][rev_f_idx] += f*1.0/num_bootstraps
										else:
												# A truly invariant site
												total_null_freq_all_snps[tp_pair_type][cohort][0] += 1.0/num_bootstraps
						
						for gene_change in gains:
							gene_name = gene_change[0]
							f = 0 if gene_name not in gene_freq_map else gene_freq_map[gene_name]
							f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()
							total_freq_gains[tp_pair_type][cohort][f_idx] += 1
						
						for gene_change in losses:
							gene_name = gene_change[0]
							f = 0 if gene_name not in gene_freq_map else gene_freq_map[gene_name]
							f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()
							total_freq_losses[tp_pair_type][cohort][f_idx] += 1
							
							# Bootstrap null model losses
							num_bootstraps = 10
							fs = choice(gene_freq_values, size=num_bootstraps, p=gene_freq_weights)
							for f in fs:
								f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
								total_null_freq_losses[tp_pair_type][cohort][f_idx] += 1.0/num_bootstraps
						
						for freq in freq_list:
							dNdS_by_prev[tp_pair_type][cohort]['N_opp'][freq] += total_non_opportunities
							dNdS_by_prev[tp_pair_type][cohort]['S_opp'][freq] += total_syn_opportunities
				
				# Pseudo counts
				nonsyn_diff_count += 1
				syn_diff_count += 1
				
				# dN/dS of each sample by timepoint pair type
				if syn_diff_count != 0:
					dNdS = (nonsyn_diff_count/float(total_non_opportunities))/(syn_diff_count/float(total_syn_opportunities))
					
					dNdS_by_tp_pair[cohort][tp_pair].append(dNdS)
	
	if False:
		pickle.dump(total_freq_snps, open(total_freq_snps_pickle_fn, 'w'))
		pickle.dump(total_freq_all_snps, open(total_freq_all_snps_pickle_fn, 'w'))
		pickle.dump(total_null_freq_snps, open(total_null_freq_snps_pickle_fn, 'w'))
		pickle.dump(total_null_freq_all_snps, open(total_null_freq_all_snps_pickle_fn, 'w'))
		pickle.dump(pooled_snp_change_distribution, open(pooled_snp_pickle_fn, 'w'))
		pickle.dump(pooled_gene_change_distribution, open(pooled_gene_pickle_fn, 'w'))
		pickle.dump(dNdS_by_prev, open(dNdS_by_prev_pickle_fn, 'w'))
		pickle.dump(pooled_snp_change_by_tp_pair, open(snp_change_by_tp_pair_fn, 'w'))
		pickle.dump(total_freq_gains, open(total_freq_gains_fn, 'w'))
		pickle.dump(total_freq_losses, open(total_freq_losses_fn, 'w'))
		pickle.dump(total_null_freq_losses, open(total_null_freq_losses_fn, 'w'))

###################################
#
# Set up figure 1: Distribution of SNV changes by variant type
# Two rows for mother-infant vs infant-infant pairs
# Three columns for the datasets
#
###################################

def is_nan(x):
	return (x is numpy.nan or x != x)

def is_inf(x):
	return (x is numpy.inf or x != x)

fig_freq_var, ax_freq_var = plt.subplots(2, 3, figsize=(8,5))
fig_freq_var.subplots_adjust(hspace=.5)

cohort_pos_map = {'backhed': (0,0), 'ferretti': (0,1), 'yassour': (0, 2)}

if sweep_type == 'full':
	cohort_ymax_map_ii = {'backhed': 80, 'ferretti': 10, 'yassour': 30}
	cohort_ymax_map_mi = {'backhed': 30, 'ferretti': 10, 'yassour': 30}
	if modification_difference_threshold == 50:
		cohort_ymax_map_ii = {'backhed': 200, 'ferretti': 50, 'yassour': 100}
		cohort_ymax_map_mi = {'backhed': 50, 'ferretti': 20, 'yassour': 50}
elif sweep_type == 'partial':
	cohort_ymax_map_ii = {'backhed': 240, 'ferretti': 50, 'yassour': 100}
	cohort_ymax_map_mi = {'backhed': 60, 'ferretti': 10, 'yassour': 40}
	if modification_difference_threshold == 50:
		cohort_ymax_map_ii = {'backhed': 300, 'ferretti': 100, 'yassour': 300}
		cohort_ymax_map_mi = {'backhed': 100, 'ferretti': 100, 'yassour': 100}
	
x = []

for cohort in ['backhed', 'yassour', 'ferretti']:
	pos1 = cohort_pos_map[cohort]
	pos2 = (pos1[0] + 1, pos1[1])
	
	for i, j in [pos1, pos2]:
		ax_freq_var[i][j].spines['top'].set_visible(False)
		ax_freq_var[i][j].spines['right'].set_visible(False)
		ax_freq_var[i][j].get_xaxis().tick_bottom()
		ax_freq_var[i][j].get_yaxis().tick_left()
		
		ax_freq_var[i][j].set_xlabel('Derived allele prevalence\nacross hosts')
		ax_freq_var[i][j].set_ylabel('# SNV changes')
		
		ax_freq_var[i][j].set_xticks(derived_virtual_xticks)
		ax_freq_var[i][j].set_xticklabels(derived_virtual_xticklabels) #,rotation='vertical')
		
		if (i, j) == pos1:
			cohort_ymax_map = cohort_ymax_map_ii
			tp_pair_type = 'II'
			if cohort == 'backhed':
				ax_freq_var[i][j].annotate("Infant-infant", xy=(0, 0.5), xytext=(-ax_freq_var[i][j].yaxis.labelpad - 5, 0), xycoords=ax_freq_var[i][j].yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
		elif (i, j) == pos2:
			cohort_ymax_map = cohort_ymax_map_mi
			tp_pair_type = 'MI'
			if cohort == 'backhed':
				ax_freq_var[i][j].annotate("Mother-infant", xy=(0, 0.5), xytext=(-ax_freq_var[i][j].yaxis.labelpad - 5, 0), xycoords=ax_freq_var[i][j].yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
		
		ax_freq_var[i][j].set_ylim([0,cohort_ymax_map[cohort]])
		ax_freq_var[i][j].bar(derived_virtual_freqs, total_freq_snps[tp_pair_type][cohort]['4D'],width=0.3,linewidth=0,facecolor='#b3de69',label='syn (4D)',zorder=3)
		
		ax_freq_var[i][j].bar(derived_virtual_freqs, total_freq_snps[tp_pair_type][cohort]['1D']+total_freq_snps[tp_pair_type][cohort]['4D'],width=0.3,linewidth=0,facecolor='#ff7f00',label='non (1D)',zorder=2)
		
		'''
		dNdS_means = {freq: [] for freq in derived_virtual_freqs}
		for f_idx in derived_virtual_freqs:
			for tp_pair in dNdS_by_tp_pair[cohort]:
				dNdS_means[f_idx].append(numpy.mean(dNdS_by_tp_pair[cohort][tp_pair][f_idx]))
		
		dNdS_mean_strs = {}
		for f_idx in derived_virtual_freqs:
			dNdS_mean_strs[f_idx] = "%.2f" % numpy.mean(dNdS_means[f_idx])
				
		for x in derived_virtual_freqs:
			for y in total_freq_snps[tp_pair_type][cohort]['1D']+total_freq_snps[tp_pair_type][cohort]['4D']:
				ax_freq_var[i][j].text(x + 0.4, y + 3, "dN/dS:\n" + dNdS_mean_strs[x])
		'''
		
		bars = ax_freq_var[i][j].bar(derived_virtual_freqs, total_freq_all_snps[tp_pair_type][cohort],width=0.3,linewidth=0,facecolor='#b15928',label='(2D & 3D)',zorder=1)
		
		ax_freq_var[i][j].bar(derived_virtual_freqs-0.3, total_null_freq_all_snps[tp_pair_type][cohort],width=0.3,linewidth=0,facecolor='0.7',label='de novo\nexpectation',zorder=0)
		
		dNdS_total = 0
		
		for k in range(len(bars)):
			rect = bars[k]
			nonsyn_diff_count = dNdS_by_prev[tp_pair_type][cohort]['N'][k]
			syn_diff_count = dNdS_by_prev[tp_pair_type][cohort]['S'][k]
			total_syn_opportunities = dNdS_by_prev[tp_pair_type][cohort]['S_opp'][k]
			total_non_opportunities = dNdS_by_prev[tp_pair_type][cohort]['N_opp'][k]
			
			dNdS = (nonsyn_diff_count/float(total_non_opportunities))/(syn_diff_count/float(total_syn_opportunities))
			
			if is_nan(dNdS) or is_inf(dNdS):
				dNdS_label = '-'
			else:
				dNdS_total += dNdS
				dNdS_label = '%0.2f' % dNdS
			
			ax_freq_var[i][j].text(rect.get_x() + rect.get_width()/2.0, rect.get_height(), dNdS_label, ha='center', va='bottom', fontsize=7, color='blue')
		
		dNdS_avg = dNdS_total/len(bars)
		x.append(dNdS_avg)
		
		ax_freq_var[i][j].legend(loc='upper right',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)
		
		ax_freq_var[i][j].set_title(cohort[0].upper() + cohort[1:], loc='center')

###################################
#
# Set up figure 2: Distribution of SNV changes by timepoint pair
# Two rows for mother-infant vs infant-infant pairs
# Three columns for the datasets
#
###################################

fig_freq_tp, ax_freq_tp = plt.subplots(2, 3, figsize=(8,5))

fig_freq_tp.subplots_adjust(hspace=.5)

cohort_pos_map = {'backhed': (0,0), 'ferretti': (0,1), 'yassour': (0,2)}

if sweep_type == 'full':
	cohort_ymax_map_ii = {'backhed': 80, 'ferretti': 10, 'yassour': 30}
	cohort_ymax_map_mi = {'backhed': 30, 'ferretti': 10, 'yassour': 30}
	if modification_difference_threshold == 50:
		cohort_ymax_map_ii = {'backhed': 200, 'ferretti': 50, 'yassour': 100}
		cohort_ymax_map_mi = {'backhed': 50, 'ferretti': 20, 'yassour': 50}
elif sweep_type == 'partial':
	cohort_ymax_map_ii = {'backhed': 240, 'ferretti': 50, 'yassour': 100}
	cohort_ymax_map_mi = {'backhed': 60, 'ferretti': 10, 'yassour': 40}
	if modification_difference_threshold == 50:
		cohort_ymax_map_ii = {'backhed': 300, 'ferretti': 100, 'yassour': 300}
		cohort_ymax_map_mi = {'backhed': 100, 'ferretti': 100, 'yassour': 100}

n_colormap = cmx.get_cmap('flag', 23)
n_colors = [n_colormap(x) for x in numpy.array([x for x in range(0,23)])/23.0]

y_colormap = cmx.get_cmap('gist_ncar', 14)
y_colors = [y_colormap(x) for x in numpy.array([x for x in range(0,14)])/14.0]

for cohort in ['backhed', 'yassour', 'ferretti']:
	pos1 = cohort_pos_map[cohort]
	pos2 = (pos1[0] + 1, pos1[1])
	
	for i, j in [pos1, pos2]:	
		ax_freq_tp[i][j].spines['top'].set_visible(False)
		ax_freq_tp[i][j].spines['right'].set_visible(False)
		ax_freq_tp[i][j].get_xaxis().tick_bottom()
		ax_freq_tp[i][j].get_yaxis().tick_left()
		
		ax_freq_tp[i][j].set_xlabel('Derived allele prevalence\nacross hosts')
		ax_freq_tp[i][j].set_ylabel('# SNV changes')
		
		ax_freq_tp[i][j].set_xticks(derived_virtual_xticks)
		ax_freq_tp[i][j].set_xticklabels(derived_virtual_xticklabels) #,rotation='vertical')
		
		colors = y_colors if cohort == 'yassour' else n_colors
		
		if (i, j) == pos1:
			cohort_ymax_map = cohort_ymax_map_ii
			tp_pair_type = 'II'
			if cohort == 'backhed':
				ax_freq_tp[i][j].annotate("Infant-infant", xy=(0, 0.5), xytext=(-ax_freq_tp[i][j].yaxis.labelpad - 5, 0), xycoords=ax_freq_tp[i][j].yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
		elif (i, j) == pos2:
			cohort_ymax_map = cohort_ymax_map_mi
			tp_pair_type = 'MI'
			if cohort == 'backhed':
				ax_freq_tp[i][j].annotate("Mother-infant", xy=(0, 0.5), xytext=(-ax_freq_tp[i][j].yaxis.labelpad - 5, 0), xycoords=ax_freq_tp[i][j].yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
		
		ax_freq_tp[i][j].set_ylim([0,cohort_ymax_map[cohort]])
		
		tp_pair_dict = {}
		for tp_pair in pooled_snp_change_by_tp_pair[cohort].keys():
			temp_type = timepoint_pair_type(tp_pair)
			if temp_type == tp_pair_type: # Exclude mom-gestation and 3mo postpartum
				tp_pair_dict[tp_pair] = pooled_snp_change_by_tp_pair[cohort][tp_pair]
		z = len(tp_pair_dict)
		prev_heights = numpy.zeros_like(derived_virtual_freqs)
		for tp_pair in tp_pair_dict.keys():	
			ax_freq_tp[i][j].bar(derived_virtual_freqs, tp_pair_dict[tp_pair] + prev_heights,width=0.3,linewidth=0,facecolor=colors[z-1],label=format_tp_label(cohort, tp_pair),zorder=z)
			prev_heights += tp_pair_dict[tp_pair]
			z -= 1
		
		ax_freq_tp[i][j].bar(derived_virtual_freqs-0.3, total_null_freq_all_snps[tp_pair_type][cohort],width=0.3,linewidth=0,facecolor='0.7',label='de novo\nexpectation',zorder=0)
		
		ax_freq_tp[i][j].legend(loc='upper right',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)
		
		ax_freq_tp[i][j].set_title(cohort[0].upper() + cohort[1:], loc='center')

###################################
#
# Set up figure 3: dN/dS by timepoint pair
#
###################################

def space_to_newline(str):
	result = ''
	for char in str:
		result += '\n' if char == ' ' else char
	return result

fig_dnds, ax_dnds = plt.subplots(1, 3, figsize=(12,3), sharey=True, gridspec_kw={'width_ratios': [1, 1.3, 2], 'wspace': 0.02})

i = 0
for cohort in ['hmp', 'backhed', 'yassour']:
	
	ax_dnds[i].set_title(cohort[0].upper() + cohort[1:])	
	data, labels = [], []
	
	for tp_pair in dNdS_by_tp_pair[cohort].keys():
		data.append(dNdS_by_tp_pair[cohort][tp_pair])
		labels.append(space_to_newline(format_tp_label(cohort, tp_pair)))
	
	j = 1	
	for subdata in data:
		mean_str = "%.2f" % (numpy.mean(subdata))
		labels[j-1] += ("\n" + 'Mean:\n' + mean_str)
		ax_dnds[i].plot(len(subdata) * [j], subdata, 'x')
		j += 1
	
	# ax_dnds[i].plot(range(len(data)), [numpy.mean(subdata) for subdata in data])
	
	ax_dnds[i].set_xticks(list(range(0, j+1)))
	ax_dnds[i].set_xticklabels([''] + labels + [''])
	ax_dnds[i].set_xlabel("Timepoint Pair")
	
	i += 1

ax_dnds[0].set_ylabel("dN/dS")

'''
dNdS_means_by_tp_pair = {cohort: {} for cohort in cohorts}
for cohort in cohorts:
	for tp_pair in dNdS_by_tp_pair[cohort]:
		dNdS_means_by_tp_pair[cohort][tp_pair] = numpy.mean(dNdS_by_tp_pair[cohort][tp_pair])
'''

###################################
#
# Set up figure 4: Gene changes
#
###################################

fig_gene_freq, ax_gene_freq = plt.subplots(2, 3, figsize=(8,5))
fig_gene_freq.subplots_adjust(hspace=.5)

cohort_pos_map = {'backhed': (0,0), 'ferretti': (0,1), 'yassour': (0, 2)}

cohort_ymax_map_ii = {'backhed': 100, 'ferretti': 10, 'yassour': 100}
cohort_ymax_map_mi = {'backhed': 30, 'ferretti': 10, 'yassour': 15}

for cohort in ['backhed', 'yassour', 'ferretti']:
	pos1 = cohort_pos_map[cohort]
	pos2 = (pos1[0] + 1, pos1[1])
	
	for i, j in [pos1, pos2]:	
		ax_gene_freq[i][j].spines['top'].set_visible(False)
		ax_gene_freq[i][j].spines['right'].set_visible(False)
		ax_gene_freq[i][j].get_xaxis().tick_bottom()
		ax_gene_freq[i][j].get_yaxis().tick_left()
		 
		ax_gene_freq[i][j].set_xlabel('Gene prevalence across hosts')
		ax_gene_freq[i][j].set_ylabel('# gene changes')

		ax_gene_freq[i][j].set_xlim([gene_freq_xticks[0],gene_freq_xticks[-1]])
		ax_gene_freq[i][j].set_xticks(gene_freq_xticks)
		ax_gene_freq[i][j].set_xticklabels(gene_freq_xticklabels) #,rotation='vertical')
		
		if (i, j) == pos1:
			cohort_ymax_map = cohort_ymax_map_ii
			tp_pair_type = 'II'
			if cohort == 'backhed':
				ax_gene_freq[i][j].annotate("Infant-infant", xy=(0, 0.5), xytext=(-ax_gene_freq[i][j].yaxis.labelpad - 5, 0), xycoords=ax_gene_freq[i][j].yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
		elif (i, j) == pos2:
			cohort_ymax_map = cohort_ymax_map_mi
			tp_pair_type = 'MI'
			if cohort == 'backhed':
				ax_gene_freq[i][j].annotate("Mother-infant", xy=(0, 0.5), xytext=(-ax_gene_freq[i][j].yaxis.labelpad - 5, 0), xycoords=ax_gene_freq[i][j].yaxis.label, textcoords='offset points', size='large', ha='right', va='center')
		
		ax_gene_freq[i][j].set_ylim([0,cohort_ymax_map[cohort]])
		
		ax_gene_freq[i][j].bar(gene_gain_virtual_freqs, total_freq_gains[tp_pair_type][cohort],width=0.3,linewidth=0,facecolor='#b3de69',label='gain')		
		ax_gene_freq[i][j].bar(gene_loss_virtual_freqs, total_freq_losses[tp_pair_type][cohort],width=0.3,linewidth=0, facecolor='#ff7f00',label='loss')		
		ax_gene_freq[i][j].bar(gene_loss_virtual_freqs-0.3, total_null_freq_losses[tp_pair_type][cohort],width=0.3,linewidth=0, facecolor='0.7',label='de novo\nexpectation')
		if (i, j) == (1, 1):
			ax_gene_freq[i][j].legend(loc='upper center',frameon=False,fontsize=5.5,numpoints=1,ncol=3,handlelength=1)
		ax_gene_freq[i][j].set_title(cohort[0].upper() + cohort[1:], loc='center')

###################################
#
# Examine genes that undergo SNP changes in multiple infants
#
###################################
'''
popular_genes = []

for gene in modified_gene_subjects.keys():
	sample_pairs = modified_gene_subjects[gene]
	if len(sample_pairs) > 1:
		popular_genes.append(gene)
		print(gene)
		print(sample_pairs)

for gene in popular_genes:
	gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
	print(gene_freq_map[gene])
'''
# Save figures

sys.stderr.write("Saving figures...\t")
fig_freq_var.savefig('%s/prevalence_latest/freq_snp_changes_by_vartype_%s_%s_%s_%s_%s.pdf' % (config.analysis_directory, prevalence_cohort, within_host_type, within_mi_pair_type, sweep_type, modification_difference_threshold),bbox_inches='tight')
fig_freq_tp.savefig('%s/prevalence_latest/freq_snp_changes_by_timepoints_%s_%s_%s_%s_%s.pdf' % (config.analysis_directory, prevalence_cohort, within_host_type, within_mi_pair_type, sweep_type, modification_difference_threshold),bbox_inches='tight')
# fig_dnds.savefig('%s/prevalence_latest/modification_dnds_by_timepoints_%s_%s_%s_%s.pdf' % (config.analysis_directory, prevalence_cohort, within_host_type, within_mi_pair_type, sweep_type),bbox_inches='tight')
fig_gene_freq.savefig('%s/prevalence_latest/freq_gene_changes_%s_%s_%s_%s_%s.pdf' % (config.analysis_directory, prevalence_cohort, within_host_type, within_mi_pair_type, sweep_type, modification_difference_threshold),bbox_inches='tight')
sys.stderr.write("Done!\n")

