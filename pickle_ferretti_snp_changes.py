from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils, core_gene_utils, gene_diversity_utils
import numpy as np
from numpy.random import choice, random as np_random, randint
import random
from collections import defaultdict
import pickle
import os, sys

# ======================================================
# Examines all consecutive timepoint pairs within hosts
# across all cohorts, and pickles SNP/gene change info
# ======================================================

# Parameters
sweep_type = 'full' # assume full for now
pp_prev_cohort = 'all'
min_coverage = 0
thresholds = {'full': (0.2, 0.8), 'partial': (0.35, 0.65)}
lower_threshold, upper_threshold = thresholds[sweep_type]

clade_divergence_threshold = 3e-02 # TODO: change to top level clade definition later

min_sample_size = 3
variant_types = ['1D','4D']
within_host_type = 'nonconsecutive' # consecutive timepoints
min_snp_change_sample_size = 5

# For partitioning SNVs according to prevalence
derived_freq_bins = np.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = np.arange(0,len(derived_freq_bins)-1)
derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
derived_virtual_xticklabels = ['0','.01','.1','.5','.9','.99','1']

# For partitioning genes into different prevalence classes
gene_freq_bins = np.array([-1,0.1,0.5,0.9,2])
gene_freq_xticks			= [-4, -3,	-2,		-1,		0,	 1,		 2,		3, 4]
gene_freq_xticklabels = ['0','0.1','0.5', '0.9','1','0.9','0.5', '0.1','0']
gene_gain_virtual_freqs = np.array([3.5,2.5,1.5,0.5])
gene_loss_virtual_freqs = np.array([-3.5,-2.5,-1.5,-0.5])

# Sample-subject-order maps
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Timepoint pair types
tp_pair_names = ['MM', 'MI', 'II', 'AA']

# Cohorts
cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'hmp']
mi_cohorts = ['backhed', 'ferretti', 'yassour', 'shao']

# Prevalence cohorts
prev_cohorts = ['all', 'hmp', 'infant', 'nonpremie', 'mother']

# Samples for each cohort
samples = {cohort: su.get_sample_names(cohort) for cohort in cohorts}
hmp_samples = su.get_sample_names('hmp')
mother_samples = su.get_sample_names('mother')
infant_samples = su.get_sample_names('infant')

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# ===================================================================
# Species SNP/gene change distributions

# species -> (sample1, sample2) -> # SNP differences
snp_changes = {species: {} for species in good_species_list}

# ===================================================================

for species_name in good_species_list[::-1]:
	
	sys.stderr.write("\nProcessing %s...\n" % species_name)
	
	# Do not restrict to QP samples!
	
	# Load temporal change map
	sys.stderr.write("Loading pre-computed temporal changes...\n")
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name, prev_cohort=pp_prev_cohort, min_coverage=min_coverage) # Default min coverage 20
	sys.stderr.write("Done!\n")
	
	# Load private SNV map
	private_snv_map = snps_utils.load_private_snv_map(species_name, prev_cohort=pp_prev_cohort)
	
	# Get all mother-infant comparisons
	cohort = 'ferretti'
	
	desired_samples = su.get_sample_names(cohort)
	
	# These indices are w.r.t. desired_samples
	same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)
	diff_subject_idxs = su.calculate_ordered_diff_subject_pairs(sample_order_map, desired_samples)
	
	# Loop over different pairs of within-host samples
	for sample_pair_idx in range(len(same_subject_idxs[0])):			 
			
			sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
			sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
			
			'''
			# Checks if among those samples from different hosts,
			# at least one of them has nonzero SNP and gene opportunities
			good_idxs = su.calculate_samples_in_different_subjects(sample_subject_map, combined_qp_samples, sample_i)
			good_idxs *= ((snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5))
			
			# FIRST FILTER
			if good_idxs.sum() < 1:
				sys.stderr.write("Not enough other-host samples!\n")
				continue
			'''
			
			# SNP temporal changes
			L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
			'''
			# SECOND FILTER
			if L<config.min_opportunities:
				sys.stderr.write("Not enough SNP opportunities (should be >=100,000)!\n")
				continue
			'''
			nerr = L*perr
			
			num_mutations = len(mutations)
			num_reversions = len(reversions)
			num_snp_changes = num_mutations + num_reversions
			
			# Gene temporal changes
			gene_L, gene_perr, gains, losses = temporal_changes_utils.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
			
			gene_nerr = gene_L*gene_perr
			num_gains = len(gains)
			num_losses = len(losses)
			num_gene_changes = num_gains+num_losses
			'''
			# THIRD FILTER
			if (perr<-0.5) or (gene_perr < -0.5):
				sys.stderr.write("Perr too high!\n")
				continue
			
			# FOURTH FILTER
			if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
				sys.stderr.write("Nerr too high!\n")
				continue # Only take things with low-ish FPR
			'''
			# ===============================================================
			# Store information
			# ===============================================================
			snp_changes[species_name][(sample_i, sample_j)] = num_snp_changes
				

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles/cov%i_prev_%s/nonconsecutive/" % (ddir, min_coverage, pp_prev_cohort)
os.system('mkdir -p %s' % pdir)

pickle.dump(snp_changes, open('%s/ferretti_all_snp_changes_%s.pkl' % (pdir, sweep_type), 'wb'))

sys.stderr.write("Done!\n")

