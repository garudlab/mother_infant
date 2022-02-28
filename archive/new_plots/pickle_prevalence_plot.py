from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils
import numpy as np
from numpy.random import choice, random, randint
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
good_species_list = parse_midas_data.parse_good_species_list()

# For partitioning SNVs according to prevalence
derived_freq_bins = np.array([-1,0,0.01,0.1,0.5,0.9,0.99,2])
derived_virtual_freqs = np.arange(0, len(derived_freq_bins) - 1)
def get_f_idx(f):
    return ((f>derived_freq_bins[:-1])*
            (f<=derived_freq_bins[1:])).argmax()

# Prevalence cohorts
prev_cohorts = ['all', 'infant', 'mother', 'hmp', 'premie', 'nonpremie']

# ===================================================================
# Pickles

# binned dN, dS num/denom for each species, cohort, tp pair
# total opportunities considered from no change and modification events
# prev cohort -> species -> cohort -> timepoint pair -> [syn diffs, syn opps, nonsyn diffs, nonsyn opps]
dNdS_distribution = {prev_cohort: {species: {cohort: {} for cohort in cohorts} for species in good_species_list} for prev_cohort in prev_cohorts}

# prevalence of sites in modifications or replacements
# prev cohort -> event type -> species -> cohort -> tp pair -> list of dicts: variant type -> list of prevalence values
event_types = ["modification", "replacement"]
prev_distribution = {prev_cohort: {event_type: {species: {cohort: defaultdict(list) for cohort in cohorts} for species in good_species_list} for event_type in event_types} for prev_cohort in prev_cohorts}

# This is the null distribution of prevalence (conditioned on total # of SNVs)
null_prev_distribution = {prev_cohort: {event_type: {species: {cohort: {} for cohort in cohorts} for species in good_species_list} for event_type in event_types} for prev_cohort in prev_cohorts}

# ===================================================================

def get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map):
	gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
	
	f1 = A1*1.0/D1
	f2 = A2*1.0/D2
	# if f1 is the frequency of the more prevalent allele at timepoint 1
	# if f2 is the frequency of the more 
	# and f is the prevalence of the less prevalent allele
	
	# reversion is when site goes back to common allele
	
	is_reversion = (f1>f2) # f2 is the prevalence of the less prevalent allele, always
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

def get_snv_prevbin_counts(snv_freq_map, private_snv_map, derived_freq_bins):
	bin_counts = np.zeros(len(derived_virtual_freqs))
	for location_tuple in snv_freq_map:
		f = snv_freq_map[location_tuple] if location_tuple not in private_snv_map else -0.5
		bin_counts[get_f_idx(f)] += 1 # One site with this f
		bin_counts[get_f_idx(1-f)] += 1 # Also consider alt allele
		# Equal weight is sketchy, but we can think of this as
		# double the actual number of opportunities
	return bin_counts

def normalize(arr):
	try:
		return arr/float(np.sum(arr)) # Assume numpy array
	except: # Expect error to be because everything is 0
		return arr

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
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name) # Default min coverage 20
	sys.stderr.write("Done!\n")
	
	# Load private SNV map
	private_snv_map = snps_utils.load_private_snv_map(species_name)
	
	# Load SNP prevalences
	snv_freq_map = {}
	for prev_cohort in prev_cohorts:
		# Each key is a location tuple; value is prevalence (f)
		snv_freq_map[prev_cohort] = snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=True)
	
	# Get proportion of "all" SNPs by SNP prevalence
	snv_prevbin_counts = get_snv_prevbin_counts(snv_freq_map[prev_cohort], private_snv_map, derived_freq_bins)
	snv_prevbin_counts_norm = normalize(snv_prevbin_counts)
	
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
				
				# set event type
				if num_snp_changes == 0:
					event_type = 'no change'
				elif num_snp_changes < 100:
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
							
							# Now draw a null prevalence from the genome
							snv_freq_keys = snv_freq_map[prev_cohort].keys()
							L = snp_opportunity_matrix[i,j]
							L_snv = len(snv_freq_map[prev_cohort]) # A slight overestimate
							snv_fraction = L_snv*1.0/L
							num_bootstraps = 10
							for bootstrap_idx in xrange(0,num_bootstraps):
							
								if random()<snv_fraction:
									# A polymorphic site
							
									random_snv_idx = randint(0,len(snv_freq_keys))
									random_snv_location = snv_freq_keys[random_snv_idx]
									f = snv_freq_map[prev_cohort][random_snv_location]
									
									rev_f = 1-f
									
									if random_snv_location in private_snv_map:
										# A private SNV. Use private bins
										f_idx = 0
										rev_f_idx = -1
									else:						
										f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()
										rev_f_idx = ((rev_f>derived_freq_bins[:-1])*(rev_f<=derived_freq_bins[1:])).argmax()		
									
									# Now add in probability weight
									if tp_pair not in null_prev_distribution[prev_cohort][event_type][species_name][cohort]:
										null_prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair] = np.zeros_like(derived_virtual_freqs)*1.0
									
									null_prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair][f_idx] += (1-f)*1.0/num_bootstraps
									null_prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair][rev_f_idx] += f*1.0/num_bootstraps
								
								else:
										# A truly invariant site
										if tp_pair not in null_prev_distribution[prev_cohort][event_type][species_name][cohort]:
											null_prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair] = np.zeros_like(derived_virtual_freqs)*1.0
										
										null_prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair][0] += 1.0/num_bootstraps
						
						prev_distribution[prev_cohort][event_type][species_name][cohort][tp_pair].append(fs_dict)
				
				# dNdS information
				if num_snp_changes < 100:				
					for prev_cohort in prev_cohorts:
						subdict = dNdS_distribution[prev_cohort][species_name][cohort]
						
						if tp_pair not in subdict:
							subdict[tp_pair] = { cat: np.zeros(len(derived_virtual_freqs)) for cat in ['dSnum', 'dSden', 'dNnum', 'dNden'] }
						
						syn_change_counts = np.zeros(len(derived_virtual_freqs))
						non_change_counts = np.zeros(len(derived_virtual_freqs))
						
						for snp_change in (mutations + reversions):
							gene_name, _, _, variant_type, A1, D1, A2, D2 = snp_change							
							f = get_sweep_prevalence(snp_change, snv_freq_map[prev_cohort], private_snv_map)
							if variant_type == '1D':
								syn_change_counts[get_f_idx(f)] += 1
							elif variant_type == '4D':
								non_change_counts[get_f_idx(f)] += 1
						
						subdict[tp_pair]['dSnum'] += syn_change_counts
						subdict[tp_pair]['dSden'] += (syn_opps * snv_prevbin_counts_norm)
						subdict[tp_pair]['dNnum'] += non_change_counts
						subdict[tp_pair]['dNden'] += (nonsyn_opps * snv_prevbin_counts_norm)				
				
				

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles" % ddir

pickle.dump(dNdS_distribution, open('%s/dNdS_prev_distribution_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(prev_distribution, open('%s/prev_distribution_v2_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(null_prev_distribution, open('%s/null_prev_distribution_%s.pkl' % (pdir, sweep_type), 'wb'))

sys.stderr.write("Done!\n")
