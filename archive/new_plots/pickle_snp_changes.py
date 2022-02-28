from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils
import numpy as np
from numpy.random import choice
from collections import defaultdict
import pickle
import sys

# ======================================================
# Examines all consecutive timepoint pairs within hosts
# across all cohorts, and pickles information about
# SNP and gene changes
# ======================================================

# Parameters
sweep_type = sys.argv[1]
thresholds = {'full': (0.2, 0.8), 'partial': (0.35, 0.65)}
lower_threshold, upper_threshold = thresholds[sweep_type]

min_sample_size = 3
variant_types = ['1D','4D']
within_host_type = 'consecutive' # consecutive timepoints
min_snp_change_sample_size = 5

# Sample-subject-order maps
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Timepoint pair types
tp_pair_names = ['MM', 'MI', 'II', 'AA']

# Cohorts
cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'olm', 'hmp']
mi_cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'olm']

# Samples for each cohort
samples = {cohort: su.get_sample_names(cohort) for cohort in cohorts}
hmp_samples = su.get_sample_names('hmp')
mother_samples = su.get_sample_names('mother')
infant_samples = su.get_sample_names('infant')

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# For partitioning SNVs according to prevalence
derived_freq_bins = np.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = np.arange(0, len(derived_freq_bins) - 1)

# ===================================================================
# Species SNP/gene change distributions

# cohort -> species -> tp_pair -> list of (# SNP changes, perr, L)
species_snp_changes = {cohort: defaultdict(dict) for cohort in cohorts}
species_between_snp_changes = {cohort: defaultdict(dict) for cohort in cohorts}

# cohort -> species -> tp_pair -> list of (# gene gains, # gene losses, perr, L)
species_gene_changes = {cohort: defaultdict(dict) for cohort in cohorts}

# cohort -> species -> tp_pair -> list of (# gene changes, perr, L)
species_between_gene_changes = {cohort: defaultdict(dict) for cohort in cohorts}

# TP pair type -> cohort -> variant type -> # SNP changes by prevalence class
total_freq_snps = {tpp_type: {cohort: \
									{vt: np.zeros_like(derived_virtual_freqs) for vt in variant_types} \
									for cohort in cohorts} for tpp_type in tp_pair_names}

# Same, but holds expectation (relative fraction of opportunities) for 1D/4D
# SNP changes, conditioned on prevalence
total_null_freq_snps =	{tpp_type: {cohort: \
												{vt: np.zeros_like(derived_virtual_freqs) for vt in variant_types} \
												for cohort in cohorts} for tpp_type in tp_pair_names} 

# syn and nonsyn diffs/opps for each species, cohort, tp pair
# event type -> species -> cohort -> timepoint pair -> [syn diffs, syn opps, nonsyn diffs, nonsyn opps
event_types = ["no change", "modification", "replacement"]
dNdS_distribution = {event_type: {species: {cohort: defaultdict(list) for cohort in cohorts} for species in good_species_list} for event_type in event_types}

# prevalence of sites in modifications or replacements
event_types = ["modification", "replacement"]
prev_cohorts = ['all', 'infant', 'mother', 'hmp', 'premie', 'nonpremie']
prev_distribution = {prev_cohort: {event_type: {species: {cohort: defaultdict(list) for cohort in cohorts} for species in good_species_list} for event_type in event_types} for prev_cohort in prev_cohorts}

# ===================================================================

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

for species_name in good_species_list[::-1]:
	
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
	
	# Load substitution rates for all QP samples
	sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
	substitution_rate_map = substitution_rates_utils.load_substitution_rate_map(species_name)
	
	if substitution_rate_map == {}: # Not enough haploid samples
		sys.stderr.write("Not enough haploid samples!\n")
		continue
	
	sys.stderr.write("Calculating SNV matrix...\n")
	dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = substitution_rates_utils.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=combined_qp_samples)
	
	snp_difference_matrix = snp_mut_difference_matrix + snp_rev_difference_matrix
	snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
	snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
	sys.stderr.write("Done!\n")
	
	sys.stderr.write("Loading gene matrix...\n")
	gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = substitution_rates_utils.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=combined_qp_samples)
	gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
	gene_opportunity_matrix = gene_loss_opportunity_matrix
	gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
	sys.stderr.write("Done!\n")
	
	sys.stderr.write("Loading 1D & 4D opportunity matrices...\n")	
	difference_matrices, opportunity_matrices = {}, {}	
	for var_type in variant_types:		
		matrix_samples, difference_matrix, opportunity_matrix = substitution_rates_utils.calculate_matrices_from_substitution_rate_map(substitution_rate_map, var_type, allowed_samples=combined_qp_samples)		
		difference_matrices[var_type] = difference_matrix
		opportunity_matrices[var_type] = opportunity_matrix
	
	sys.stderr.write("Done!\n")
	
	# Load temporal change map
	sys.stderr.write("Loading pre-computed temporal changes...\n")
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name)
	sys.stderr.write("Done!\n")
	
	# Load private SNV map
	private_snv_map = snps_utils.load_private_snv_map(species_name)
	
	# Load SNP prevalences
	snv_freq_map = {}
	for prev_cohort in prev_cohorts:
		snv_freq_map[prev_cohort] = snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=True)
	
	# Loop over different cohorts
	for cohort in cohorts:		
		desired_samples = qp_sample_lists[cohort]
		
		same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type, one_per_mi_pair=False)
		
		# Loop over different pairs of within-host samples
		for sample_pair_idx in range(len(same_subject_idxs[0])):			 
				
				sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
				sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
				tp_pair = su.sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples)
				
				i = combined_sample_idx_map[sample_i]
				j = combined_sample_idx_map[sample_j]
				
				# Checks if among those samples from different hosts,
				# at least one of them has nonzero SNP and gene opportunities
				good_idxs = su.calculate_samples_in_different_subjects(sample_subject_map, combined_qp_samples, sample_i)
				good_idxs *= ( (snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5) )
				
				if good_idxs.sum() < 1:
					sys.stderr.write("Not enough other-host samples!\n")
					continue
				
				matrix_idx_i = matrix_samples.index(sample_i)
				matrix_idx_j = matrix_samples.index(sample_j)
				
				# Numbers of site differences and opportunities between the timepoints
				nonsyn_diffs = difference_matrices['1D'][matrix_idx_i][matrix_idx_j]
				nonsyn_opps = opportunity_matrices['1D'][matrix_idx_i][matrix_idx_j]		
				syn_diffs = difference_matrices['4D'][matrix_idx_i][matrix_idx_j]
				syn_opps = opportunity_matrices['4D'][matrix_idx_i][matrix_idx_j]
				
				# SNP temporal changes
				L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
				
				if L<config.min_opportunities:
					sys.stderr.write("Not enough SNP opportunities (should be >=100,000)!\n")
					continue
				
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
				
				if (perr<-0.5) or (gene_perr < -0.5):
					sys.stderr.write("Perr too high!\n")
					continue
				
				if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
					sys.stderr.write("Nerr too high!\n")
					continue # Only take things with low-ish FPR
				
				# Store information
				
				snp_change_tup = (num_snp_changes, perr, L)
				between_snp_change_tup = (choice(snp_difference_matrix[i, good_idxs]), perr, L)
				
				try:
					species_snp_changes[cohort][species_name][tp_pair].append(snp_change_tup)
				except:
					species_snp_changes[cohort][species_name][tp_pair] = []
					species_snp_changes[cohort][species_name][tp_pair].append(snp_change_tup)
				
				try:
					species_between_snp_changes[cohort][species_name][tp_pair].append(between_snp_change_tup)
				except:
					species_between_snp_changes[cohort][species_name][tp_pair] = []
					species_between_snp_changes[cohort][species_name][tp_pair].append(between_snp_change_tup)
				
				gene_change_tup = (num_gains, num_losses, gene_perr, gene_L)
				between_gene_change_tup = (choice(gene_difference_matrix[i, good_idxs]), gene_perr, gene_L)
				
				try:
					species_gene_changes[cohort][species_name][tp_pair].append(gene_change_tup)
				except:
					species_gene_changes[cohort][species_name][tp_pair] = []
					species_gene_changes[cohort][species_name][tp_pair].append(gene_change_tup)
				
				try:
					species_between_gene_changes[cohort][species_name][tp_pair].append(between_gene_change_tup)
				except:
					species_between_gene_changes[cohort][species_name][tp_pair] = []
					species_between_gene_changes[cohort][species_name][tp_pair].append(between_gene_change_tup)
				
				# set event type
				if num_snp_changes == 0:
					event_type = 'no change'
				elif num_snp_changes < 20:
					event_type = 'modification'
				elif num_snp_changes > 400:
					event_type = 'replacement'
				else:
					continue
				
				# prevalence information
				if num_snp_changes > 0:
					for prev_cohort in prev_cohorts:
						
						fs_dict = defaultdict(list)
						
						for snp_change in (mutations + reversions):
							gene_name, _, _, variant_type, A1, D1, A2, D2 = snp_change							
							f = get_sweep_prevalence(snp_change, snv_freq_map[prev_cohort], private_snv_map)
							fs_dict[variant_type].append(f)
						
						prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair].append(fs_dict)
				
				# dNdS information
				subdict = dNdS_distribution[event_type][species_name][cohort]
				if tp_pair not in subdict:
					subdict[tp_pair] = [nonsyn_diffs, nonsyn_opps, syn_diffs, syn_opps]
				else:
					subdict[tp_pair][0] += nonsyn_diffs
					subdict[tp_pair][1] += nonsyn_opps
					subdict[tp_pair][2] += syn_diffs
					subdict[tp_pair][3] += syn_opps

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles" % ddir

pickle.dump(dNdS_distribution, open('%s/dNdS_distribution_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(species_snp_changes, open('%s/species_snp_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(species_between_snp_changes, open('%s/species_between_snp_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(species_gene_changes, open('%s/species_gene_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(species_between_gene_changes, open('%s/species_between_gene_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(prev_distribution, open('%s/prev_distribution_%s.pkl' % (pdir, sweep_type), 'wb'))

sys.stderr.write("Done!\n")
