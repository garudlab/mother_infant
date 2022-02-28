from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils, midas_db_utils, core_gene_utils, parse_patric
import numpy as np
from numpy.random import choice, random, randint
from collections import defaultdict
import pickle
import sys

# ======================================================
# Examines all consecutive timepoint pairs within hosts
# across all cohorts, and pickles information about
# SNP and gene changes.
# Takes in one argument: full or partial
# ======================================================

# Parameters
sweep_type = sys.argv[1]
thresholds = {'full': (0.2, 0.8), 'partial': (0.35, 0.65)}
lower_threshold, upper_threshold = thresholds[sweep_type]

min_sample_size = 3 # at least this many QP samples to consider species
variant_types = ['1D','4D'] # what sites to include in syn/opp matrices
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

# Modification SNP changes
# species > cohort > subject > contig-position tuple >
# list of snp change tuples of the new form
# (tp_pair, gene_name, variant_type, A1, D1, A2, D2)
snp_modifications_by_site = {species: {cohort: {} for cohort in cohorts} for species in good_species_list}

# Modification genes
# species > cohort > subject > sample[tp] pair > gene_name > count
snp_modification_genes = {species: {cohort: {} for cohort in cohorts} for species in good_species_list}

# Modification gene gains
# species > cohort > subject > sample[tp] pair > gene_name > count
modification_gene_gains = {species: {cohort: {} for cohort in cohorts} for species in good_species_list}

# Modification gene losses
# species > cohort > subject > sample[tp] pair > gene_name > count
modification_gene_losses = {species: {cohort: {} for cohort in cohorts} for species in good_species_list}

# Modification gene gains - present genes null
# species > cohort > subject > sample[tp] pair > gene_name > count
modification_gene_gains_present_null = {species: {cohort: {} for cohort in cohorts} for species in good_species_list}

# Modification gene losses - present genes null
# species > cohort > subject > sample[tp] pair > gene_name > count
modification_gene_losses_present_null = {species: {cohort: {} for cohort in cohorts} for species in good_species_list}

# ===================================================================

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
	
	# Load data for nulls
	# Note: gene_samples has same content as combined_qp_samples,
	# but order may vary...
	sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
	gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name, allowed_samples=combined_qp_samples)
	sys.stderr.write("Done!\n")
	
	# Later, consider gene present if its copy number 0.5-2
	# Matrix rows - genes, cols - gene_samples
	gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
	
	# Load temporal change map
	sys.stderr.write("Loading pre-computed temporal changes...\n")
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name, 20) # Default min coverage 20
	sys.stderr.write("Done!\n")
	
	# Load private SNV map
	private_snv_map = snps_utils.load_private_snv_map(species_name)
	
	# Get gene information
	genome_ids = midas_db_utils.get_ref_genome_ids(species_name)
	non_shared_genes = core_gene_utils.parse_non_shared_pangenome_genes(species_name)
	gene_desc = parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
	
	# Loop over different cohorts
	for cohort in cohorts:		
		qp_samples_cohort = qp_sample_lists[cohort]
		
		# Note: only consecutive infant timepoints are included,
		# but all mother-infant timepoints combos are included
		same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, qp_samples_cohort, within_host_type=within_host_type, one_per_mi_pair=False)
		
		# Loop over different pairs of within-host samples
		for sample_pair_idx in range(len(same_subject_idxs[0])):			 
				
				sample_i = qp_samples_cohort[same_subject_idxs[0][sample_pair_idx]] 
				sample_j = qp_samples_cohort[same_subject_idxs[1][sample_pair_idx]]
				# print(sample_order_map[sample_i][0] + " " + str(sample_order_map[sample_i][1]) + " > " + sample_order_map[sample_j][0] + " " + str(sample_order_map[sample_j][1]))
				tp_pair = su.sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples)
				
				i = combined_sample_idx_map[sample_i]
				j = combined_sample_idx_map[sample_j]
				
				subject = sample_subject_map[sample_i]
				if subject != sample_subject_map[sample_j]:
					print("Subject weirdness for samples " + sample_i + ", " + sample_j)
				
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
				
				# Modification event genes and reversion information
				if event_type == 'modification':
					
					# SNP changes
					for snp_change in (mutations + reversions):						
						gene_id, contig, position, variant_type, A1, D1, A2, D2 = snp_change
						gene_name = gene_desc[gene_id] if gene_id in gene_desc else 'N/A'
						
						if subject not in snp_modification_genes[species_name][cohort]:
							snp_modification_genes[species_name][cohort][subject] = {}
							snp_modification_genes[species_name][cohort][subject][tp_pair] = defaultdict(int)
						elif tp_pair not in snp_modification_genes[species_name][cohort][subject]:
							snp_modification_genes[species_name][cohort][subject][tp_pair] = defaultdict(int)
						
						if subject not in snp_modifications_by_site[species_name][cohort]:
							snp_modifications_by_site[species_name][cohort][subject] = defaultdict(list)
						
						snp_modification_genes[species_name][cohort][subject][tp_pair][(gene_id, gene_name)] += 1
						snp_modifications_by_site[species_name][cohort][subject][(contig, position)].append((tp_pair, (gene_id, gene_name), variant_type, A1, D1, A2, D2))
					
					# Gene gains
					for gene_gain in gains:
						gene_id, D1, Dm1, D2, Dm2 = gene_gain
						gene_name = gene_desc[gene_id] if gene_id in gene_desc else 'N/A'
						
						if subject not in modification_gene_gains[species_name][cohort]:
							modification_gene_gains[species_name][cohort][subject] = {}
						
						if tp_pair not in modification_gene_gains[species_name][cohort][subject]:
							modification_gene_gains[species_name][cohort][subject][tp_pair] = defaultdict(int)
						
						modification_gene_gains[species_name][cohort][subject][tp_pair][(gene_id, gene_name)] += 1
					
					# Gene losses
					for gene_loss in losses:
						gene_id, D1, Dm1, D2, Dm2 = gene_loss
						gene_name = gene_desc[gene_id] if gene_id in gene_desc else 'N/A'
						
						if subject not in modification_gene_losses[species_name][cohort]:
							modification_gene_losses[species_name][cohort][subject] = {}
						
						if tp_pair not in modification_gene_losses[species_name][cohort][subject]:
							modification_gene_losses[species_name][cohort][subject][tp_pair] = defaultdict(int)
						
						modification_gene_losses[species_name][cohort][subject][tp_pair][(gene_id, gene_name)] += 1
					
					# Null: pick same number (as total number of changes) of genes
					# randomly from all genes present at either time point
					gene_samples_idx1 = list(gene_samples).index(sample_i)
					gene_samples_idx2 = list(gene_samples).index(sample_j)

					present_gene_idxs = []
					for idx in [gene_samples_idx1, gene_samples_idx2]:
						present_gene_idxs.extend(np.nonzero((gene_copynum_matrix[:,idx]>0.5)*(gene_copynum_matrix[:,idx]<2))[0])
					
					num_trials = 50
					
					for trial in range(num_trials):
						try:						
							present_gene_gain_null_idxs = choice(present_gene_idxs, num_gains)
							present_gene_loss_null_idxs = choice(present_gene_idxs, num_losses)
							
							for idx in present_gene_gain_null_idxs:
								gene_id = gene_names[idx]
								gene_name = gene_desc[gene_id] if gene_id in gene_desc else 'N/A'
																
								if subject not in modification_gene_gains_present_null[species_name][cohort]:
									modification_gene_gains_present_null[species_name][cohort][subject] = {}
								
								if tp_pair not in modification_gene_gains_present_null[species_name][cohort][subject]:
									modification_gene_gains_present_null[species_name][cohort][subject][tp_pair] = defaultdict(int)
								
								modification_gene_gains_present_null[species_name][cohort][subject][tp_pair][(gene_id, gene_name)] += (1.0/num_trials)
							
							for idx in present_gene_loss_null_idxs:
								gene_id = gene_names[idx]
								gene_name = gene_desc[gene_id] if gene_id in gene_desc else 'N/A'
																
								if subject not in modification_gene_losses_present_null[species_name][cohort]:
									modification_gene_losses_present_null[species_name][cohort][subject] = {}
								
								if tp_pair not in modification_gene_losses_present_null[species_name][cohort][subject]:
									modification_gene_losses_present_null[species_name][cohort][subject][tp_pair] = defaultdict(int)
								
								modification_gene_losses_present_null[species_name][cohort][subject][tp_pair][(gene_id, gene_name)] += (1.0/num_trials)				
						except:
							print('something went wrong for')
							print(species_name)
							print(subject)
							continue
					

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles" % ddir

pickle.dump(snp_modification_genes, open('%s/snp_modification_genes_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(snp_modifications_by_site, open('%s/snp_modifications_by_site_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(modification_gene_gains, open('%s/modification_gene_gains_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(modification_gene_gains_present_null, open('%s/modification_gene_gains_present_null_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(modification_gene_losses, open('%s/modification_gene_losses_%s.pkl' % (pdir, sweep_type), 'wb'))
pickle.dump(modification_gene_losses_present_null, open('%s/modification_gene_losses_present_null_%s.pkl' % (pdir, sweep_type), 'wb'))

sys.stderr.write("Done!\n")
