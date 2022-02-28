# Investigating gene gain/loss vs alpha diversity
# For HMP adults only

import matplotlib	 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sample_utils as su, config, parse_midas_data, stats_utils, sfs_utils
import pylab, sys, numpy as np, random, math
import calculate_temporal_changes
from collections import defaultdict
import bz2

def list_to_colors(input_list):
	cmap = matplotlib.cm.hsv
	
	input_list_dict = {}
	i = 0
	for elem in set(input_list):
		input_list_dict[elem] = i
		i += 1
	
	norm = matplotlib.colors.Normalize(vmin=0, vmax=i-1)
	m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
	
	color_list = []
	for elem in input_list:
		color_i = input_list_dict[elem]
		color_list.append(m.to_rgba(color_i))
	
	return color_list

def colors_to_legend_elements(colors, labels):
	legend_elements = {}
	for color, label in zip(colors, labels):
		legend_elements[color] = matplotlib.patches.Patch(facecolor=color,                         label=label)
	return legend_elements.values()

#######################################################
# Standard header to read in argument information
#######################################################
import argparse
parser = argparse.ArgumentParser()
# parser.add_argument('--species', type=str, help='Run the script for one specified species')
parser.add_argument('--sweep-type', type=str, help="Full or partial sweep", default="partial")

args = parser.parse_args()
sweep_type = args.sweep_type

if sweep_type not in ['full', 'partial']:
	sys.exit("Invalid sweep-type. Choose from full, partial")

if sweep_type == 'full':
	lower_threshold, upper_threshold = 0.2, 0.8
elif sweep_type == 'partial':
	lower_threshold, upper_threshold = 0.35, 0.65
#######################################################

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.calculate_sample_subject_map(subject_sample_map)
sys.stderr.write("Done!\n")

# Samples
all_samples = sample_country_map.keys() 
hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

# Species list
good_species_list = parse_midas_data.parse_good_species_list()

# Relative abundance file
relab_fpath = "%s/species/relative_abundance.txt.bz2" % (config.data_directory)
relab_file = open(relab_fpath, 'r')
decompressor = bz2.BZ2Decompressor()
raw = decompressor.decompress(relab_file.read())
data = [row.split('\t') for row in raw.split('\n')]
data.pop() # Get rid of extra element due to terminal newline
header = data[0]

# Dictionary: host -> gene gain/loss numbers aggregated across species
host_gene_gain_dict = defaultdict(int)
host_gene_loss_dict = defaultdict(int)

# Dictionary: host-species pair -> gene gain/loss numbers
host_species_gene_gain_dict = defaultdict(int)
host_species_gene_loss_dict = defaultdict(int)

# Dictionary: host -> alpha diversity at (first, second) timepoint
# Only include if modification event happened in host (longest)
host_alpha_diversity_dict = {}

# Generate alpha diversity dictionary
alpha_div_dict = {}
for i in range(1, len(header)):
	acc = 0
	for row in data[1:]:
		rel_ab = float(row[i])
		if rel_ab != 0:
			acc += (rel_ab * math.log(rel_ab))
	alpha_div_dict[header[i]] = (acc*-1)

# Missing samples in relative abundance file?
bad_samples = []

# Modification / replacement event thresholds
modification_difference_threshold = config.modification_difference_threshold
replacement_difference_threshold = config.replacement_difference_threshold

# Store temporal change info
for species_name in good_species_list:
	
	print("Working on " + species_name)
		
	desired_samples = sorted(su.calculate_qp_samples(hmp_samples, species_name)['qp'])
	_, same_subject_idxs, _ = su.calculate_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type='longest')	
	temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
	
	for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
		
		sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
		sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
		
		subject = sample_subject_map[sample_i]
		
		gene_opps, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
		all_changes = gains + losses
		
		# Compute number of gene gains/losses
		num_gains = len(gains) if gains != None else 0
		num_losses = len(losses) if losses != None else 0
		num_gene_changes = num_gains + num_losses
		
		# Don't want to look at things with high error rates
		if (gene_perr<-0.5) or (gene_perr>0.5):
			continue
		
		# Store gene change information
		host_gene_gain_dict[subject] += num_gains
		host_gene_loss_dict[subject] += num_losses
		host_species_gene_gain_dict[(subject, species_name)] = num_gains
		host_species_gene_loss_dict[(subject, species_name)] = num_losses
		
		# Store alpha (Shannon) diversity for earlier sample
		if sample_i in alpha_div_dict and sample_j in alpha_div_dict:
			host_alpha_diversity_dict[subject] = (alpha_div_dict[sample_i], alpha_div_dict[sample_j])
		elif sample_i in alpha_div_dict: # Only first timepoint has alpha div info
			bad_samples.append(sample_i)
			host_alpha_diversity_dict[subject] = (alpha_div_dict[sample_i], None)
		elif sample_j in alpha_div_dict: # Only second timepoint has alpha div info
			bad_samples.append(sample_j)
			host_alpha_diversity_dict[subject] = (None, alpha_div_dict[sample_j])
		else: # Neither available
			bad_samples.append(sample_i)
			bad_samples.append(sample_j)

# Check out which samples have no relative abundance info
print(bad_samples)

# Store information in .csv
gene_changes_csv = open('%s/hmp_gene_changes.csv' % (config.analysis_directory), 'w')

gene_changes_csv.write('subject_id,species_id,num_gene_gains,num_gene_losses,alpha_div_tp1,alpha_div_tp2\n')

for host, species in host_species_gene_gain_dict:
	gene_gains = host_species_gene_gain_dict[(host, species)]
	gene_losses = host_species_gene_loss_dict[(host, species)]
	
	try:
		alpha_div1, alpha_div2 = host_alpha_diversity_dict[host]
		alpha_div1 = alpha_div1 if alpha_div1 != None else 'NA'
		alpha_div2 = alpha_div2 if alpha_div2 != None else 'NA'
	except:
		continue
	
	gene_changes_csv.write('%s,%s,%s,%s,%s,%s\n' % (host, species, gene_gains, gene_losses, alpha_div1, alpha_div2))

gene_changes_csv.close()

# Set up gene gain/loss data for plot
host_alpha_diversities = []
host_gene_gain_counts = []
host_gene_loss_counts = []

host_alpha_div_diffs = []
host_gene_gain_counts_for_diffs = []
host_gene_loss_counts_for_diffs = []

for subject in host_alpha_diversity_dict:
	alpha_div1, alpha_div2 = host_alpha_diversity_dict[subject]
	gain_count = host_gene_gain_dict[subject]
	loss_count = host_gene_loss_dict[subject]
	
	if alpha_div1 != None:
		host_alpha_diversities.append(alpha_div1)
		host_gene_gain_counts.append(gain_count)
		host_gene_loss_counts.append(loss_count)
	
	if alpha_div1 != None and alpha_div2 != None:
		alpha_div_diff = alpha_div2 - alpha_div1
		host_alpha_div_diffs.append(alpha_div_diff)
		host_gene_gain_counts_for_diffs.append(gain_count)
		host_gene_loss_counts_for_diffs.append(loss_count)

host_gene_change_counts = np.array(host_gene_gain_counts) + np.array(host_gene_loss_counts)
host_gene_change_counts_for_diffs = np.array(host_gene_gain_counts_for_diffs) + np.array(host_gene_loss_counts_for_diffs)

host_species_gene_gain_alpha_divs = []
host_species_gene_gain_counts = []
host_species_gene_gain_species = []

host_species_gene_gain_alpha_div_diffs = []
host_species_gene_gain_counts_for_diffs = []
host_species_gene_gain_species_for_diffs = []

host_species_gene_loss_alpha_divs = []
host_species_gene_loss_counts = []
host_species_gene_loss_species = []

host_species_gene_loss_alpha_div_diffs = []
host_species_gene_loss_counts_for_diffs = []
host_species_gene_loss_species_for_diffs = []

for subject, species in host_species_gene_gain_dict:
	
	gain_count = host_species_gene_gain_dict[(subject, species)]
	
	if subject in host_alpha_diversity_dict:
		alpha_div1, alpha_div2 = host_alpha_diversity_dict[subject]
		
		if alpha_div1 != None:
			host_species_gene_gain_alpha_divs.append(alpha_div1)
			host_species_gene_gain_counts.append(gain_count)
			host_species_gene_gain_species.append(species)
		
		if alpha_div1 != None and alpha_div2 != None:
			alpha_div_diff = alpha_div2 - alpha_div1
			host_species_gene_gain_alpha_div_diffs.append(alpha_div_diff)
			host_species_gene_gain_counts_for_diffs.append(gain_count)
			host_species_gene_gain_species_for_diffs.append(species)

for subject, species in host_species_gene_loss_dict:
	
	loss_count = host_species_gene_loss_dict[(subject, species)]
	
	if subject in host_alpha_diversity_dict:
		alpha_div1, alpha_div2 = host_alpha_diversity_dict[subject]
		
		if alpha_div1 != None:
			host_species_gene_loss_alpha_divs.append(alpha_div1)
			host_species_gene_loss_counts.append(loss_count)
			host_species_gene_loss_species.append(species)
		
		if alpha_div1 != None and alpha_div2 != None:
			alpha_div_diff = alpha_div2 - alpha_div1
			host_species_gene_loss_alpha_div_diffs.append(alpha_div_diff)
			host_species_gene_loss_counts_for_diffs.append(loss_count)
			host_species_gene_loss_species_for_diffs.append(species)

# Plot 1: one point per host, aggregate across species
fig, ax = plt.subplots()
ax.plot(host_alpha_diversities, host_gene_gain_counts, 'ro', label='Gene gains')
# ax.set_yscale('log')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene gains, aggregated across species")
ax.set_title("Alpha diversity vs number of gene gains, HMP adults (%i subjects)" % len(host_gene_gain_counts))

fig.savefig('%s/alpha_div_vs_gene_gains_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_diversities, host_gene_loss_counts, 'bo', label='Gene losses')
# ax.set_yscale('log')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene losses, aggregated across species")
ax.set_title("Alpha diversity vs number of gene losses, HMP adults (%i subjects)" % len(host_gene_loss_counts))
fig.savefig('%s/alpha_div_vs_gene_losses_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_diversities, host_gene_change_counts, 'go', label='Combined')
# ax.set_yscale('log')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene changes (gain or loss), aggregated across species")
ax.set_title("Alpha diversity vs number of gene changes, HMP adults (%i subjects)" % len(host_gene_change_counts))

fig.savefig('%s/alpha_div_vs_gene_changes_hmp_agg.png' % (config.analysis_directory))

# =================================================================
# Plot 2: one point per host-species pair
# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_gain_species)
ax.scatter(host_species_gene_gain_alpha_divs, host_species_gene_gain_counts, c=colors, edgecolors='none')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene gains per host-species pair")
ax.set_title("Alpha div. vs # gene gains,\nHMP adults (n=%i)" % len(host_species_gene_gain_counts))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_gain_species), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_vs_gene_gains_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')

# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_loss_species)
ax.scatter(host_species_gene_loss_alpha_divs, host_species_gene_loss_counts, c=colors, edgecolors='none')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of gene losses per host-species pair")
ax.set_title("Alpha div. vs # gene losses,\nHMP adults (n=%i)" % len(host_species_gene_loss_counts))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_loss_species), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_vs_gene_loss_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')

# =================================================================
# Plot 3: one point per host, alpha div diffs
# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_div_diffs, host_gene_gain_counts_for_diffs, 'ro', label='Gene gains')
# ax.set_yscale('log')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene gains, aggregated across species")
ax.set_title("Alpha diversity change vs # gene gains,\nHMP adults (%i subjects)" % len(host_gene_gain_counts_for_diffs))

plt.grid(True,which="both", linestyle='--')
fig.savefig('%s/alpha_div_diff_vs_gene_gains_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_div_diffs, host_gene_loss_counts_for_diffs, 'bo', label='Gene losses')
# ax.set_yscale('log')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene losses, aggregated across species")
ax.set_title("Alpha diversity change vs # gene losses,\nHMP adults (%i subjects)" % len(host_gene_loss_counts_for_diffs))

fig.savefig('%s/alpha_div_diff_vs_gene_losses_hmp_agg.png' % (config.analysis_directory))

# =================================================================

fig, ax = plt.subplots()
ax.plot(host_alpha_div_diffs, host_gene_change_counts_for_diffs, 'go', label='Combined')
# ax.set_yscale('log')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene changes (gain or loss), aggregated across species")
ax.set_title("Alpha diversity change vs # gene changes,\nHMP adults (%i subjects)" % len(host_gene_change_counts_for_diffs))

fig.savefig('%s/alpha_div_diff_vs_gene_changes_hmp_agg.png' % (config.analysis_directory))

# =================================================================
# Plot 4: one point per host-species pair, alpha div diffs
# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_gain_species_for_diffs)
ax.scatter(host_species_gene_gain_alpha_div_diffs, host_species_gene_gain_counts_for_diffs, c=colors, edgecolors='none')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene gains per host-species pair")
ax.set_title("Alpha div. vs # gene gains,\nHMP adults (%i host-species pairs)" % len(host_species_gene_gain_counts_for_diffs))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_gain_species_for_diffs), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_diff_vs_gene_gains_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')

# =================================================================

fig, ax = plt.subplots(figsize=(8,12))
colors = list_to_colors(host_species_gene_loss_species_for_diffs)
ax.scatter(host_species_gene_loss_alpha_div_diffs, host_species_gene_loss_counts_for_diffs, c=colors, edgecolors='none')
ax.set_xlabel("Change in Shannon diversity")
ax.set_ylabel("Number of gene losses per host-species pair")
ax.set_title("Alpha div. vs # gene losses,\nHMP adults (%i host-species pairs)" % len(host_species_gene_loss_counts_for_diffs))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_gene_loss_species_for_diffs), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_diff_vs_gene_loss_hmp_ind.png' % (config.analysis_directory), bbox_inches='tight')