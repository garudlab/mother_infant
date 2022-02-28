# Investigating 'amount of evolution' vs alpha diversity
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

# Dictionary: host -> number of SNP changes aggregated across species
# Only include if modification event happened in host (longest)
host_snp_change_dict = defaultdict(int)

# Dictionary: host-species pair -> number of SNP changes
host_species_snp_change_dict = defaultdict(int)

# Dictionary: host -> alpha diversity at first timepoint
# Only include if modification event happened in host (longest)
host_alpha_diversity_dict = {}

# Dictionary: none / mod / replace -> host -> number of events [i.e. species with the relevant event type, for longest timepoint pair]
host_change_type_dict = {type: defaultdict(int) for type in ['none', 'mod', 'replace']}

# Dictionary: sample(-species pair) -> within-sample polymorphism (average estimate)
sample_species_polymorphism_dict = defaultdict(dict)

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

# Modification event threshold
modification_difference_threshold = config.modification_difference_threshold

# Replacement event threshold
replacement_difference_threshold = config.replacement_difference_threshold

# Store temporal change info
for species_name in good_species_list:
	
	print("Working on " + species_name)
	
	# Store within-host polymorphism for ALL samples
	_, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
	for sample in hmp_samples:		
		if sample not in sfs_map: # TODO
			continue
		within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
		try:
			within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05) # 95% CI
		except:
			continue
		sample_species_polymorphism_dict[sample][species_name] = (within_rate_lower + within_rate_upper)/2.0
	
	desired_samples = sorted(su.calculate_qp_samples(hmp_samples, species_name)['qp'])
	_, same_subject_idxs, _ = su.calculate_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type='longest')	
	temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)	
	
	for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
		
		sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
		sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
		
		L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold, upper_threshold)
		
		nerr = L*perr
		num_mutations = len(mutations)
		num_reversions = len(reversions)
		num_snp_changes = num_mutations+num_reversions
		
		subject = sample_subject_map[sample_i]
		
		# Ignore iffy sample pairs
		if L < config.min_opportunities:
			continue
		if (perr<-0.5) or (perr>0.5):
			continue					
		if (nerr > max([0.5, 0.1*num_snp_changes])):
			continue # Only take things with low-ish FPR
		if num_snp_changes == 0:
			host_change_type_dict['none'][subject] += 1
		elif num_snp_changes >= replacement_difference_threshold:
			host_change_type_dict['replace'][subject] += 1
		else:
			# Count anything not replacement and nonzero as modification here
			host_change_type_dict['mod'][subject] += 1
		
		if num_snp_changes >= modification_difference_threshold:
			continue
		
		# Store number of SNP changes in modification event for host
		host_snp_change_dict[subject] += num_snp_changes
		host_species_snp_change_dict[(subject, species_name)] = num_snp_changes
		
		# Store alpha (Shannon) diversity for earlier sample
		if sample_i in alpha_div_dict:
			host_alpha_diversity_dict[subject] = alpha_div_dict[sample_i]
		elif sample_j in alpha_div_dict:
			bad_samples.append(sample_i)
			host_alpha_diversity_dict[subject] = alpha_div_dict[sample_j]
		else:
			bad_samples.append(sample_j)

# Check out which samples have no relative abundance info
print(bad_samples)

# Set up proportion replacement vs modification for plot
# Threshold on at least 10 species present
prop_replace = []
prop_mod = []
prop_repmod_alpha_diversities = []

for subject in host_alpha_diversity_dict:
	
	if subject not in host_change_type_dict['none']:
		none_count = 0
	else:
		none_count = host_change_type_dict['none'][subject]
	
	if subject not in host_change_type_dict['replace']:
		replace_count = 0
	else:
		replace_count = host_change_type_dict['replace'][subject]
	
	if subject not in host_change_type_dict['mod']:
		mod_count = 0
	else:
		mod_count = host_change_type_dict['mod'][subject]
	
	total_count = replace_count + mod_count + none_count
	
	if total_count > 5:
		prop_replace.append(float(replace_count)/total_count)
		prop_mod.append(float(mod_count)/total_count)
		prop_repmod_alpha_diversities.append(host_alpha_diversity_dict[subject])

# Set up modification SNP change data for plot
host_alpha_diversities = []
host_snp_change_counts = []

for subject in host_alpha_diversity_dict:
	host_alpha_diversities.append(host_alpha_diversity_dict[subject])
	host_snp_change_counts.append(host_snp_change_dict[subject])

host_species_alpha_diversities = []
host_species_snp_change_counts = []
host_species_species = []

for subject, species in host_species_snp_change_dict:
	if subject in host_alpha_diversity_dict:
		host_species_alpha_diversities.append(host_alpha_diversity_dict[subject])
		host_species_snp_change_counts.append(host_species_snp_change_dict[(subject, species)])
		host_species_species.append(species)

# Set up within-sample polymorphism data for plot
sample_species_alpha_diversities = []
sample_species_polymorphisms = []

sample_alpha_diversities = []
sample_polymorphisms = []

for sample in sample_species_polymorphism_dict:
	if sample in alpha_div_dict:
		sample_alpha_diversities.append(alpha_div_dict[sample])
		agg_polymorphism = 0
		for species in sample_species_polymorphism_dict[sample]:
			sample_species_alpha_diversities.append(alpha_div_dict[sample])
			sample_species_polymorphisms.append(sample_species_polymorphism_dict[sample][species])
			agg_polymorphism += sample_species_polymorphism_dict[sample][species]
		sample_polymorphisms.append(agg_polymorphism/len(sample_species_polymorphism_dict[sample]))

# Plot 0: sample-species alpha diversity vs. polymorphism
fig, ax = plt.subplots()
ax.scatter(sample_species_alpha_diversities, sample_species_polymorphisms)
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Within-sample polymorphism per sample-species pair")
ax.set_title("Alpha diversity vs polymorphism, HMP adults (n=%i)" % len(sample_species_polymorphisms))

fig.savefig('%s/alpha_div_vs_polymorphism_hmp.png' % (config.analysis_directory), bbox_inches='tight')

# Plot 0.5: sample alpha diversity vs. polymorphism (averaged across species)
fig, ax = plt.subplots()
ax.scatter(sample_alpha_diversities, sample_polymorphisms)
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Within-sample polymorphism per sample, averaged across species")
ax.set_title("Alpha diversity vs polymorphism, HMP adults (n=%i)" % len(sample_species_polymorphisms))

fig.savefig('%s/alpha_div_vs_polymorphism_avg_hmp.png' % (config.analysis_directory), bbox_inches='tight')

# Plot 1: one point per host, aggregate across species
fig, ax = plt.subplots()
ax.scatter(host_alpha_diversities, host_snp_change_counts)
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of SNP changes, aggregated across species")
ax.set_title("Alpha diversity vs number of SNP changes, HMP adults (n=%i)" % len(host_snp_change_counts))

fig.savefig('%s/alpha_div_vs_snp_changes_hmp_%s_agg.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')

# Plot 2: one point per host-species pair
fig, ax = plt.subplots()
colors = list_to_colors(host_species_species)
ax.scatter(host_species_alpha_diversities, host_species_snp_change_counts, c=colors, edgecolors='none')
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Number of SNP changes per host-species pair")
ax.set_title("Alpha div. vs # SNP changes, HMP adults (n=%i)" % len(host_species_snp_change_counts))

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', handles=colors_to_legend_elements(colors, host_species_species), fontsize='x-small', bbox_to_anchor=(1, 0.5))

fig.savefig('%s/alpha_div_vs_snp_changes_hmp_%s_ind.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')

# Plot 3: Alpha diversity vs. fraction modification
fig, ax = plt.subplots()
ax.scatter(prop_repmod_alpha_diversities, prop_mod, color='r', label="Modification")
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Fraction of present species which are modified")
ax.set_title("Alpha diversity vs. fraction of modification events, HMP adults (n=%i)" % len(prop_mod))

fig.savefig('%s/alpha_div_vs_prop_mod_hmp_%s.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')

# Plot 4: Alpha diversity vs. fraction replacement
fig, ax = plt.subplots()
ax.scatter(prop_repmod_alpha_diversities, prop_replace, color='b', label="Replacement")
ax.set_xlabel("Shannon diversity")
ax.set_ylabel("Fraction of present species which are replaced")
ax.set_title("Alpha diversity vs. fraction of replacement events, HMP adults (n=%i)" % len(prop_replace))

fig.savefig('%s/alpha_div_vs_prop_replace_hmp_%s.png' % (config.analysis_directory, sweep_type), bbox_inches='tight')
