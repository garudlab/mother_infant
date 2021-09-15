from utils import config, parse_midas_data, sample_utils as su, temporal_changes_utils, stats_utils, midas_db_utils, parse_patric, snps_utils
from collections import defaultdict
import math, random, numpy as np
import pickle, sys, bz2, os
import matplotlib.pyplot as plt

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

# Sample-subject-order maps
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Cohorts
cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'hmp']
mi_cohorts = ['backhed', 'ferretti', 'yassour', 'shao']

# Samples for each cohort
samples = {cohort: su.get_sample_names(cohort) for cohort in cohorts}
hmp_samples = su.get_sample_names('hmp')
mother_samples = su.get_sample_names('mother')
infant_samples = su.get_sample_names('infant')
olm_samples = su.get_sample_names('olm')
infant_samples_no_olm = [sample for sample in infant_samples if sample not in olm_samples]

# Sample-timepoint map
mi_sample_day_dict = su.get_mi_sample_day_dict(exclude_cohorts=['olm'])

# Consider all mother/infant samples here
desired_samples = mother_samples + infant_samples
# These indices are w.r.t. desired_samples
same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type='consecutive', one_per_mi_pair=False, infant_timepoint_pref='first')
diff_subject_idxs = su.calculate_ordered_diff_subject_pairs(sample_order_map, desired_samples)

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# Relative abundance file
relab_fpath = "%s/species/relative_abundance.txt.bz2" % (config.data_directory)
relab_file = open(relab_fpath, 'r')
decompressor = bz2.BZ2Decompressor()
raw = decompressor.decompress(relab_file.read())
data = [row.split('\t') for row in raw.split('\n')]
data.pop() # Get rid of extra element due to terminal newline
header = su.parse_merged_sample_names(data[0]) # species_id, samples...

# Load species presence/absence information
sample_species_dict = defaultdict(set)

for row in data[1:]:
    species = row[0]
    for relab_str, sample in zip(row[1:], header[1:]):
        relab = float(relab_str)
        if relab > 0:
            sample_species_dict[sample].add(species)

# Custom mother-infant test
is_mi = lambda sample_i, sample_j: ((sample_i in mother_samples and sample_j in infant_samples_no_olm) and mi_sample_day_dict[sample_i] >= 0 and mi_sample_day_dict[sample_j] <= 7)

# ===================================================================
# (sample1, sample2) -> species -> (shared/not_shared, # SNP differences)
# if not_shared, # SNP differences should be -1
transmission_info = defaultdict(dict)
# ===================================================================

for species_name in good_species_list[::-1]:
	
	sys.stderr.write("\nProcessing %s...\n" % species_name)
	
	# Do not restrict to QP samples! However, keep high coverage filter later
	
	# Load temporal change map
	sys.stderr.write("Loading pre-computed temporal changes...\n")
	temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name, prev_cohort=pp_prev_cohort, min_coverage=min_coverage) # Default min coverage 20
	sys.stderr.write("Done!\n")
	
	# Loop over different pairs of within-host samples
	for sample_pair_idx in range(len(same_subject_idxs[0])):			 
			
			sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
			sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
			
			if not is_mi(sample_i, sample_j):
				continue
			
			# Shared species
			shared_species = sample_species_dict[sample_i].intersection(sample_species_dict[sample_j])
			
			# SNP temporal changes
			L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
			
			# FILTER
			if L<config.min_opportunities:
				sys.stderr.write("Not enough SNP opportunities (should be >=100,000)!\n")
				continue
			
			nerr = L*perr			
			num_mutations = len(mutations)
			num_reversions = len(reversions)
			num_snp_changes = num_mutations + num_reversions
			
			# ===============================================================
			# Store information
			# ===============================================================
			shared_status = 'shared' if species_name in shared_species else 'not_shared'
			transmission_info[(sample_i, sample_j)][species_name] = (shared_status, num_snp_changes)
				

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles/cov%i_prev_%s/" % (ddir, min_coverage, pp_prev_cohort)
os.system('mkdir -p %s' % pdir)

pickle.dump(transmission_info, open('%s/transmission_info.pkl' % (pdir), 'wb'))

sys.stderr.write("Done!\n")

# Get dem numbers

num_transmissions = 0
num_total = 0
num_shared_species_all_vals = []
num_shared_species_highcov_vals = []
shared_species_highcov_agg = []
half_prop_transmissions = []

for sample_i, sample_j in transmission_info:
	
	# Shared species
	shared_species = sample_species_dict[sample_i].intersection(sample_species_dict[sample_j])
	
	num_shared_species_all = len(shared_species)
	num_shared_species_highcov = 0
	num_transmission = 0
	num_replacement = 0
	shared_species_highcov = set()
	
	for species in transmission_info[(sample_i, sample_j)]:
		
		shared_status, num_snp_changes = transmission_info[(sample_i, sample_j)][species]
		
		if shared_status == 'shared':
			num_shared_species_highcov += 1
		
		shared_species_highcov.add(species)
		num_total += 1
		
		if num_snp_changes <= 20: # Transmission
			num_transmission += 1
			num_transmissions += 1
		elif num_snp_changes > 20: # Replacement
			num_replacement += 1
	
	shared_species_highcov_agg += list(shared_species_highcov)
	prop_transmission = float(num_transmission)/float(num_shared_species_highcov)
	half_prop_transmission = float(num_transmission)/float(2*num_shared_species_highcov)
	half_prop_transmissions.append(half_prop_transmission)
	num_shared_species_all_vals.append(num_shared_species_all)
	num_shared_species_highcov_vals.append(num_shared_species_highcov)

shared_species_highcov_count = defaultdict(int)

for species in shared_species_highcov_agg:
	shared_species_highcov_count[species] += 1
