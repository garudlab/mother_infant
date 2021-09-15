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

# ===================================================================
# Species SNP/gene change distributions

# species -> (sample1, sample2) -> [list of SNP change tuples]
# OR # SNP changes if that number exceeds 20 (replacement)
# where each SNP change tuple consists of
# (gene_name, contig, position, variant_type, A1, D1, A2, D2)
snp_changes = {species: {} for species in good_species_list}

# species -> (sample1, sample2) -> [nonsyn_diffs, nonsyn_opps, syn_diffs, syn_opps]
dnds_info = {species: {} for species in good_species_list}

# species -> (sample1, sample2) -> (gene gain tuples, gene loss tuples)
# OR (# gene gains, # gene losses) if replacement event
# where each gene change tuple consists of...
gene_changes = {species: {} for species in good_species_list}

# species -> (sample1, sample2) -> [list of (variant_type, prevalence, opps) tuples]
# only for modification events
snp_change_freqs = {species: defaultdict(list) for species in good_species_list}
# species -> (sample1, sample2) -> prev_cohort -> [list of (prevalence, weight) tuples]
snp_change_null_freqs = {species: defaultdict(dict) for species in good_species_list}

# species -> (sample1, sample2) -> [list of prevalences]
# only for modification events
gene_gain_freqs = {species: defaultdict(list) for species in good_species_list}
gene_loss_freqs = {species: defaultdict(list) for species in good_species_list}
gene_loss_null_freqs = {species: defaultdict(list) for species in good_species_list}

# species -> (sample1, sample2) -> random number of SNP/gene differences between sample1 and unrelated host
# note that sample2 doesn't really mean anything here
between_snp_change_counts = {species: {} for species in good_species_list}
between_gene_change_counts = {species: {} for species in good_species_list}

# 3 null distributions of: for each SNP change...
# present gene: assign random gene that is present at either timepoint
# between host: assign random gene that changes between hosts
# pangenome: assign random gene from pangenome
# ...
# species -> (sample1, sample2) -> gene ID -> # occurrence (normalized by num_trials)
snp_change_present_gene_null = {species: defaultdict(dict) for species in good_species_list}
snp_change_between_host_null = {species: defaultdict(dict) for species in good_species_list}
snp_change_pangenome_null = {species: defaultdict(dict) for species in good_species_list}

# ===================================================================

for species_name in good_species_list[::-1]:
	
	sys.stderr.write("\nProcessing %s...\n" % species_name)
	
	# Grab QP samples for this species
	qp_sample_lists = {}
	for cohort in cohorts:
		qp_sample_lists[cohort] = sorted(su.load_qp_samples(samples[cohort], species_name, prev_cohort=pp_prev_cohort)['qp'])
	
	combined_qp_samples = sorted(su.flatten([qp_sample_lists[cohort] for cohort in cohorts]))
	combined_sample_name_idx_map = {combined_qp_samples[i] : i for i in range(len(combined_qp_samples))}
	
	# Using all QP samples to threshold on sample size
	if len(combined_qp_samples) < min_sample_size:
		sys.stderr.write("Not enough haploid samples!\n")
		continue
	
	# Load substitution rates for all QP samples
	sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
	substitution_rate_map = substitution_rates_utils.load_substitution_rate_map(species_name, prev_cohort=pp_prev_cohort)
	
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
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name, prev_cohort=pp_prev_cohort, min_coverage=min_coverage) # Default min coverage 20
	sys.stderr.write("Done!\n")
	
	# Load private SNV map
	private_snv_map = snps_utils.load_private_snv_map(species_name, prev_cohort=pp_prev_cohort)
	
	# Load prevalences
	snv_freq_map = {prev_cohort: snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=True) for prev_cohort in prev_cohorts}
	
	# Load gene frequencies
	gene_freq_map = {}
	for prev_cohort in prev_cohorts:
		gene_freqs = core_gene_utils.parse_gene_freqs(species_name, prev_cohort=prev_cohort)
		if len(gene_freqs) != 0:
			gene_freq_map[prev_cohort] = gene_freqs
	gene_freq_values = {prev_cohort: np.array(gene_freq_map[prev_cohort].values()) for prev_cohort in gene_freq_map}
	gene_freq_weights = {prev_cohort: (gene_freq_values[prev_cohort]*1.0/gene_freq_values[prev_cohort].sum()) for prev_cohort in gene_freq_map}
	
	# Load info for pangenome null
	pangenome_gene_names, pangenome_new_species_names = parse_midas_data.load_pangenome_genes(species_name)
	
	# Load data for between-host changes null
	gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=combined_qp_samples)
	gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
	
	# Loop over different cohorts
	for cohort in cohorts:		
		desired_samples = qp_sample_lists[cohort]
		
		# These indices are w.r.t. desired_samples
		same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type, one_per_mi_pair=False)
		diff_subject_idxs = su.calculate_ordered_diff_subject_pairs(sample_order_map, desired_samples)
		
		# These indices are w.r.t. combined_qp_samples, so that we can
		# properly index into opportunity matrices
		combined_sample_idx_map = su.calculate_sample_idx_map(desired_samples, combined_qp_samples)
		combined_same_subject_idxs = su.apply_sample_index_map_to_indices(combined_sample_idx_map, same_subject_idxs)
		combined_diff_subject_idxs = su.apply_sample_index_map_to_indices(combined_sample_idx_map, diff_subject_idxs)
		
		# These indices are w.r.t. gene_samples, so that we can
		# properly index into pangenome data
		gene_sample_idx_map = su.calculate_sample_idx_map(desired_samples, gene_samples)
		gene_same_subject_idxs = su.apply_sample_index_map_to_indices(gene_sample_idx_map, same_subject_idxs)
		gene_diff_subject_idxs = su.apply_sample_index_map_to_indices(gene_sample_idx_map, diff_subject_idxs)
		
		# This specifically gets all genes that differ between arbitrary hosts
		# within this cohort, as long as they each have a QP sample with sufficient
		# marker coverage and low enough SNP substitution rate
		between_host_gene_idxs = []
		for idx in range(len(combined_diff_subject_idxs[0])):
			combined_i = combined_diff_subject_idxs[0][idx]
			combined_j = combined_diff_subject_idxs[1][idx]
			gene_i = gene_diff_subject_idxs[0][idx]
			gene_j = gene_diff_subject_idxs[1][idx]
			if (marker_coverages[gene_i]>min_coverage) and (marker_coverages[gene_j]>min_coverage):
				if snp_substitution_rate[combined_i, combined_j] < clade_divergence_threshold:
					gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(combined_i,combined_j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
					between_host_gene_idxs.extend(gene_idxs)
		
		# Loop over different pairs of within-host samples
		for sample_pair_idx in range(len(same_subject_idxs[0])):			 
				
				sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
				sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
						
				sample_i_combined_idx = combined_same_subject_idxs[0][sample_pair_idx]
				sample_j_combined_idx = combined_same_subject_idxs[1][sample_pair_idx]
				
				sample_i_gene_idx = gene_same_subject_idxs[0][sample_pair_idx]
				sample_j_gene_idx = gene_same_subject_idxs[1][sample_pair_idx]
				
				subject = (sample_subject_map[sample_i], sample_subject_map[sample_j])
				tp_pair = su.sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples)
				
				i = combined_sample_name_idx_map[sample_i]
				j = combined_sample_name_idx_map[sample_j]
				
				matrix_idx_i = matrix_samples.index(sample_i)
				matrix_idx_j = matrix_samples.index(sample_j)
				
				# Set up data for present gene null
				# Obtain the indices of all genes which have copynum 0.5-2.0
				# at either timepoint (duplicate if show up in both samples)
				present_gene_idxs = []
				present_gene_idxs.extend(np.nonzero((gene_copynum_matrix[:,sample_i_gene_idx]>0.5)*(gene_copynum_matrix[:,sample_i_gene_idx]<2))[0])
				present_gene_idxs.extend(np.nonzero((gene_copynum_matrix[:,sample_j_gene_idx]>0.5)*(gene_copynum_matrix[:,sample_j_gene_idx]<2))[0])
				
				# Checks if among those samples from different hosts,
				# at least one of them has nonzero SNP and gene opportunities
				good_idxs = su.calculate_samples_in_different_subjects(sample_subject_map, combined_qp_samples, sample_i)
				good_idxs *= ((snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5))
				
				# FIRST FILTER
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
				
				# SECOND FILTER
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
				
				# THIRD FILTER
				if (perr<-0.5) or (gene_perr < -0.5):
					sys.stderr.write("Perr too high!\n")
					continue
				
				# FOURTH FILTER
				if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
					sys.stderr.write("Nerr too high!\n")
					continue # Only take things with low-ish FPR
				
				# ===============================================================
				# Store information
				# ===============================================================
				'''
				# dNdS information
				dnds_info[species_name][(sample_i, sample_j)] = [nonsyn_diffs, nonsyn_opps, syn_diffs, syn_opps]
				
				# Nulls
				num_trials = 100
				normalized_count = (1/float(num_trials))
				snp_change_present_gene_null[species_name][(sample_i, sample_j)] = defaultdict(int)
				snp_change_between_host_null[species_name][(sample_i, sample_j)] = defaultdict(int)
				snp_change_pangenome_null[species_name][(sample_i, sample_j)] = defaultdict(int)
				
				# If no genes available for a null, skip this QP pair
				if len(between_host_gene_idxs) > 0 and len(present_gene_idxs) > 0:				
					
					for trial in range(0,num_trials):
						# In each trial, get a random list of genes that change according to respective null
						present_gene_null_names = gene_names[choice(present_gene_idxs, num_snp_changes)]
						between_gene_null_names = gene_names[choice(between_host_gene_idxs, num_snp_changes)]					
						pangenome_gene_null_names = choice(list(pangenome_gene_names), num_snp_changes)
						
						# Then add counts normalized by total number of trails
						for gene in present_gene_null_names:
							snp_change_present_gene_null[species_name][(sample_i, sample_j)][gene] += normalized_count
						
						for gene in between_gene_null_names:
							snp_change_between_host_null[species_name][(sample_i, sample_j)][gene] += normalized_count
						
						for gene in pangenome_gene_null_names:
							snp_change_pangenome_null[species_name][(sample_i, sample_j)][gene] += normalized_count
				
				# Number of SNP differences with random unrelated host
				between_snp_change_counts[species_name][(sample_i, sample_j)] = choice(snp_difference_matrix[i, good_idxs])
				between_gene_change_counts[species_name][(sample_i, sample_j)] = choice(gene_difference_matrix[i, good_idxs])
				'''
				if num_snp_changes <= 20: # Modification event
					
					snp_changes[species_name][(sample_i, sample_j)] = (mutations + reversions) #!
					gene_changes[species_name][(sample_i, sample_j)] = (gains, losses) #!
					
					# Delve into prevalence of modified SNPs and genes
					null_freq_dict = defaultdict(list)
					
					for snp_change in (mutations + reversions):
						
						variant_type = snp_change[3]
						
						# Construct freq_dict (prevalence of sweeping allele)
						freq_dict = {}
						for prev_cohort in prev_cohorts:
							f = get_sweep_prevalence(snp_change, snv_freq_map[prev_cohort], private_snv_map)
							freq_dict[prev_cohort] = f
						
						# Construct opp_dict (number of 1D/4D opportunities for the QP pair)
						opp_dict = {'1D': opportunity_matrices['1D'][i,j], '4D': opportunity_matrices['4D'][i,j]}
						
						# Store prevalence and opportunity info
						snp_change_freqs[species_name][(sample_i, sample_j)].append((variant_type, freq_dict, opp_dict)) #!
						
						# Draw null prevalence from genome
						for prev_cohort in prev_cohorts:
							
							snv_freq_keys = snv_freq_map[prev_cohort].keys()
							
							L = snp_opportunity_matrix[i,j]
							L_snv = len(snv_freq_map[prev_cohort]) # A slight overestimate
							snv_fraction = L_snv*1.0/L
							num_bootstraps = 10
							
							for _ in range(num_bootstraps):
								
								if np_random() < snv_fraction: # Polymorphic site
									
									random_snv_location = snv_freq_keys[randint(0, len(snv_freq_keys))]
									f = 0 if random_snv_location in private_snv_map else snv_freq_map[prev_cohort][random_snv_location]
									rev_f = 1-f
									
									# Now add in probability weight
									null_freq_dict[prev_cohort].append((f, (1-f)*1.0/num_bootstraps))
									null_freq_dict[prev_cohort].append((1-f, f*1.0/num_bootstraps))
								
								else: # A truly invariant site
									
									null_freq_dict[prev_cohort].append((0, 1.0/num_bootstraps))
						
					# Store the dictionary of null prevalences by prevalence cohort
					snp_change_null_freqs[species_name][(sample_i, sample_j)] = null_freq_dict
					
					for gene_change in gains:
						
						gene_name = gene_change[0]						
						freq_dict = {}
						for prev_cohort in gene_freq_map:
							f = gene_freq_map[prev_cohort][gene_name] if gene_name in gene_freq_map[prev_cohort] else 0
							freq_dict[prev_cohort] = f
						
						gene_gain_freqs[species_name][(sample_i, sample_j)].append(freq_dict) #!
					
					for gene_change in losses:
						
						gene_name = gene_change[0]
						
						freq_dict = {}
						null_freq_dict = {}
						num_bootstraps = 10
						
						for prev_cohort in gene_freq_map:
							f = gene_freq_map[prev_cohort][gene_name] if gene_name in gene_freq_map[prev_cohort] else 0
							freq_dict[prev_cohort] = f
							
							fs = choice(gene_freq_values[prev_cohort], size=num_bootstraps, p=gene_freq_weights[prev_cohort])
							null_freq_dict[prev_cohort] = fs							
						
						gene_loss_freqs[species_name][(sample_i, sample_j)].append(freq_dict) #!
						
						gene_loss_null_freqs[species_name][(sample_i, sample_j)].append(null_freq_dict)
									
				else: # Likely replacement and too many SNPs to store info for
					snp_changes[species_name][(sample_i, sample_j)] = num_snp_changes
					gene_changes[species_name][(sample_i, sample_j)] = (num_gains, num_losses)
				

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles/cov%i_prev_%s/nonconsecutive/" % (ddir, min_coverage, pp_prev_cohort)
os.system('mkdir -p %s' % pdir)

'''
pickle.dump(snp_changes, open('%s/big_snp_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(gene_changes, open('%s/big_gene_changes_%s.pkl' % (pdir, sweep_type), 'wb'))
'''
pickle.dump(snp_change_freqs, open('%s/snp_change_freqs_with_opps_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(snp_change_null_freqs, open('%s/snp_change_null_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(gene_gain_freqs, open('%s/gene_gain_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(gene_loss_freqs, open('%s/gene_loss_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(gene_loss_null_freqs, open('%s/gene_loss_null_freqs_%s.pkl' % (pdir, sweep_type), 'wb'))
'''
pickle.dump(between_snp_change_counts, open('%s/between_snp_change_counts_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(between_gene_change_counts, open('%s/between_gene_change_counts_%s.pkl' % (pdir, sweep_type), 'wb'))

pickle.dump(snp_change_present_gene_null, open('%s/snp_change_present_gene_null.pkl' % pdir, 'wb'))
pickle.dump(snp_change_between_host_null, open('%s/snp_change_between_host_null.pkl' % pdir, 'wb'))
pickle.dump(snp_change_pangenome_null, open('%s/snp_change_pangenome_null.pkl' % pdir, 'wb'))
pickle.dump(dnds_info, open('%s/dnds_info.pkl' % pdir, 'wb'))
'''

sys.stderr.write("Done!\n")

