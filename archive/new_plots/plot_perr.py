from utils import config, sample_utils as su, parse_midas_data, sfs_utils, diversity_utils, calculate_temporal_changes
import sys, numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt

# Load subject and sample metadata
sample_order_map = su.parse_sample_order_map()

# Load good species list
good_species_list = parse_midas_data.parse_good_species_list()

# Dictionary: species -> list of perr
# Restrict to QP samples and within-host changes
all_perr_full_species = defaultdict(list)
all_perr_partial_species = defaultdict(list)

for species in good_species_list:
	
	# Load SFS map (include all variant types)
	sys.stderr.write("Loading SFS for %s...\n" % species)
	samples, sfs_map  = parse_midas_data.parse_within_sample_sfs(species)
	
	# Get QP samples
	qp_samples = list(su.calculate_qp_samples(samples, species)['qp'])
	
	final_samples = [] # hack for now
	for sample in qp_samples:
		if sample in sample_order_map:
			final_samples.append(sample)
	
	# Load temporal change map
	sys.stderr.write("Loading temporal changes for %s...\n" % species)
	temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species)
	
	# Get all within-host timepoint pairs
	same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, final_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)
	
	for i, j in zip(same_subject_idxs[0], same_subject_idxs[1]):
		
		sample_i, sample_j = final_samples[i], final_samples[j]
		
		perr_full, perr_partial = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j, dfs=[0.6, 0.3])
		
		all_perr_full_species[species].append(perr_full)
		all_perr_partial_species[species].append(perr_partial)

# Pickle!!
import pickle
pickle.dump(all_perr_full_species, open("%s/all_perr_full_species.pkl" % config.analysis_directory, 'wb'))
pickle.dump(all_perr_partial_species, open("%s/all_perr_partial_species.pkl" % config.analysis_directory, 'wb'))

# Simple histograms of full vs. partial perrs
all_perr_full = []
for species in all_perr_full_species:
	all_perr_full += all_perr_full_species[species]

all_perr_partial = []
for species in all_perr_partial_species:
	all_perr_partial += all_perr_partial_species[species]

fig, ax = plt.subplots()
ax.set_yscale('log')
ax.boxplot([all_perr_full, all_perr_partial])
ax.set_ylabel("Perr")
ax.set_xticklabels(["Full (0.6)", "Partial (0.3)"])
fig.savefig("%s/full_vs_partial_sweep_perrs.png" % config.analysis_directory)