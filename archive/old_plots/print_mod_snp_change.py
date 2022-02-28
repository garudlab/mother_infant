# Printing out all modification SNP change data

import sample_utils as su, config, parse_midas_data
import pylab, sys, numpy, random
import diversity_utils, gene_diversity_utils, sfs_utils, calculate_substitution_rates, calculate_temporal_changes, parse_patric, species_phylogeny_utils, core_gene_utils, stats_utils
from collections import defaultdict
from numpy.random import choice
import pickle

#######################################################
# Standard header to read in argument information
#######################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--species', type=str, help='Run the script for one specified species')
args = parser.parse_args()
species = args.species

if species is None:
	good_species_list = parse_midas_data.parse_good_species_list()
	pickle_fname = "%s/pickles/all_mod_snp_change_details.pkl" % config.data_directory
else:
	good_species_list = [species]
	pickle_fname = "%s/pickles/all_mod_snp_change_details/sweep0-35to0-65/all_mod_snp_change_details_%s.pkl" % (config.data_directory, species)

modification_difference_threshold = config.modification_difference_threshold # May change
min_coverage = config.min_median_coverage
clade_divergence_threshold = 1e-02 # TODO: change to top level clade definition later
num_trials=100
min_sample_size = 5

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()

hmp_samples = su.get_sample_names('hmp','all')
all_mi_samples = su.get_sample_names('all','all')
mother_samples = su.get_sample_names('mother','all')
infant_samples = su.get_sample_names('infant','all')
backhed_samples = su.get_sample_names('backhed', 'all')
ferretti_samples = su.get_sample_names('ferretti', 'all')
yassour_samples = su.get_sample_names('yassour', 'all')
sys.stderr.write("Done!\n")

# Samples to work with
cohort = hmp_samples + infant_samples + mother_samples

all_data = defaultdict(list)

for species_name in good_species_list:
	# Only plot samples above a certain depth threshold that are "haploids"
	haploid_samples = diversity_utils.calculate_haploid_samples(species_name)
	# Isolate those samples also belonging to specified set
	haploid_samples = [s for s in haploid_samples if s in cohort]
	
	temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
	
	if len(haploid_samples) < min_sample_size:
		continue
	
	same_subject_idxs = su.calculate_mi_ordered_subject_pairs(sample_order_map, haploid_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)
	
	for i, j in zip(same_subject_idxs[0], same_subject_idxs[1]):
		
		sample_i, sample_j = haploid_samples[i], haploid_samples[j]
		subject = sample_order_map[sample_i][0][:-2]
		
		try:
			snp_opportunities, snp_perr, snp_changes = temporal_change_map[(sample_i, sample_j)]['snps']
		except:
			snp_perr = -1
			print("Error: %s, %s not found in temporal change map for %s" % (sample_i, sample_j, species_name))
		
		snp_opportunities, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, 0.35, 0.65)
		snp_changes = mutations + reversions
		num_snp_changes = len(snp_changes)
		
		if num_snp_changes >= 0 and num_snp_changes <modification_difference_threshold:
			for snp_change in snp_changes:
				gene_id, contig, pos, variant_type, A1, D1, A2, D2 = snp_change
				mod_snp_change = (sample_i, sample_j, gene_id, contig, pos, variant_type, A1, D1, A2, D2, snp_perr)
				try:
					gene_name = parse_patric.load_patric_gene_name(gene_id)
					all_data[gene_name].append(mod_snp_change)
				except:
					all_data["bad_gene_id"].append(mod_snp_change)
	
pickle.dump(all_data, open(pickle_fname, "wb"))
