# This uses the postprocessed data files
import matplotlib
matplotlib.use('Agg') 
import parse_midas_data
import sample_utils
import sys
import numpy
import diversity_utils
import stats_utils
from math import log10,ceil
import matplotlib.pyplot as plt
import config
import sfs_utils

# Standard header to read in argument information
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--minfreq", help="Minimum frequency of dominant allele for polymorphism", default=0.8)
parser.add_argument("--species", help="Species to plot polymorphism for", default='Bacteroides_vulgatus_57955')
parser.add_argument("--coverage", help="Minimum coverage", default=config.min_median_coverage)

args = parser.parse_args()

species_name = args.species
species_name_clean = ' '.join(species_name.split('_')[:2])
upper_threshold = float(args.minfreq)
lower_threshold = 1 - upper_threshold
min_coverage = args.coverage

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sys.stderr.write("Done!\n") 

# Load SNP information for species_name
sys.stderr.write("Loading SFSs for %s...\t" % species_name)
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
sys.stderr.write("Done!\n")

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
samples = numpy.array(samples)

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold
desired_samples = samples[(median_coverages>=min_coverage)]
desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])

###################################
#
# Calculate within polymorphism rates
#
###################################

sample_names = []

between_rates = []

within_rates = []
within_rate_lowers = []
within_rate_uppers = []

median_depths = []
depth_lowers = []
depth_uppers = []

for sample in desired_samples:
	within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample], lower_threshold, upper_threshold)
	
	within_rate = within_sites*1.0/total_sites
	between_rate = between_sites*1.0/total_sites
	between_rates.append(between_rate)
	within_rates.append(within_rate)
	
	# Calculate 95% confidence intervals
	within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05)
	within_rate_lowers.append(within_rate_lower)
	within_rate_uppers.append(within_rate_upper)
	
	depths, counts = sfs_utils.calculate_depth_distribution_from_sfs_map(sfs_map[sample])
	dlower, dupper = stats_utils.calculate_IQR_from_distribution(depths, counts)
	dmedian = stats_utils.calculate_median_from_distribution(depths,counts)
	
	depth_lowers.append(dlower)
	depth_uppers.append(dupper)
	median_depths.append(dmedian)
	sample_names.append(sample)

# Sort them all in descending order of within-host diversity		
within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths, depth_lowers, depth_uppers, sample_names = (numpy.array(x) for x in zip(*sorted(zip(within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths,depth_lowers, depth_uppers, sample_names),reverse=True)))

within_rate_lowers = numpy.clip(within_rate_lowers, 1e-09,1)

###################################
#
# Partition data by cohorts across datasets
#
###################################

backhed_ref_freq_frac = [{} for i in range(4)]
backhed_sample_names = []
backhed_cohort_names = ['B','4M','12M','M']
for short in backhed_cohort_names:
	backhed_sample_names.append(sample_utils.get_sample_names('backhed', short))

ferretti_ref_freq_frac = [{} for i in range(6)]
ferretti_sample_names = []
ferretti_cohort_names = ['M0','I1','I2','I3','I4','I5']
for short in ferretti_cohort_names:
	ferretti_sample_names.append(sample_utils.get_sample_names('ferretti', short))

yassour_ref_freq_frac = [{} for i in range(8)]
yassour_sample_names = []
yassour_cohort_names = ['MGest', 'MBirth', 'M3', 'CBirth', 'C14', 'C1', 'C2', 'C3']
for short in yassour_cohort_names:
	yassour_sample_names.append(sample_utils.get_sample_names('yassour', short))

for rank_idx in xrange(0,len(within_rates)):
	sample = sample_names[rank_idx]
	# Check if Backhed
	for i in range(len(backhed_ref_freq_frac)):
		if sample in backhed_sample_names[i]:
			backhed_ref_freq_frac[i][sample] = within_rates[rank_idx]
	# Check if Ferretti
	for i in range(len(ferretti_sample_names)):
		if sample in ferretti_sample_names[i]:
			ferretti_ref_freq_frac[i][sample] = within_rates[rank_idx]
	# Check if Yassour
	for i in range(len(yassour_sample_names)):
		if sample in yassour_sample_names[i]:
			yassour_ref_freq_frac[i][sample] = within_rates[rank_idx]

###################################
#
# Set up figure (6 panels, 2 rows / 3 cols)
#
###################################

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharey='row',figsize=(13,6))

###################################
#
# Boxplots
#
###################################

ax1.set_ylabel("Within-sample polymorphism\n(%s $\leq f \leq$ %s)" % (lower_threshold, upper_threshold))

# Backhed
plot_data = [dict.values() for dict in backhed_ref_freq_frac]
plot_names = ["Newborns","4 Months", "12 Months", "Mothers"]
plot_names_n = ["{}\nn={}".format(a_, b_) for a_, b_ in zip(plot_names, [len(x) for x in plot_data])]

ax1.boxplot(plot_data)
ax1.set_yscale('log')
ax1.set_ylim((1e-7, 1))
ax1.tick_params(axis='x', which='major', labelsize=10)
ax1.set_xticklabels(plot_names_n)
ax1.set_title("Backhed")

# Ferretti
plot_data = [dict.values() for dict in ferretti_ref_freq_frac]
plot_names = ["Mothers","Day1", "Day3", "Week1", "Month1", "Month4"]
plot_names_n = ["{}\nn={}".format(a_, b_) for a_, b_ in zip(plot_names, [len(x) for x in plot_data])]

ax2.boxplot(plot_data)
ax2.set_yscale('log')
ax2.set_ylim((1e-7, 1))
ax2.tick_params(axis='x', which='major', labelsize=10)
ax2.set_xticklabels(plot_names_n)
ax2.set_title("%s\nFerretti" % species_name_clean)

# Yassour
plot_data = [dict.values() for dict in yassour_ref_freq_frac]
plot_names = ["M:gest","M:del","M:mon3","C:birth","C:wk2","C:mon1","C:mon2","C:mon3"]
plot_names_n = ["{}\nn={}".format(a_, b_) for a_, b_ in zip(plot_names, [len(x) for x in plot_data])]

ax3.boxplot(plot_data)
ax3.set_yscale('log')
ax3.set_ylim((1e-7, 1))
ax3.tick_params(axis='x', which='major', labelsize=9)
ax3.set_xticklabels(plot_names_n)
ax3.set_title("Yassour")

###################################
#
# What about for individuals/same subject?
#
###################################

subject_sample_map = sample_utils.parse_subject_sample_map()
sample_subject_map = sample_utils.calculate_sample_subject_map(subject_sample_map)

# Backhed
backhed_subject_rffrac_dict = {}
for i in range(len(backhed_ref_freq_frac)):
	short = backhed_cohort_names[i]
	for sample in backhed_ref_freq_frac[i].keys():
		subject = sample_subject_map[sample]
		if subject not in backhed_subject_rffrac_dict:
			backhed_subject_rffrac_dict[subject] = {}
		backhed_subject_rffrac_dict[subject][short] = backhed_ref_freq_frac[i][sample]

# 0: mother-birth, 1: birth-month4, 2: month4-month12
backhed_rffrac_diff = [[] for i in range(3)]

for subject in backhed_subject_rffrac_dict:
	subdict = backhed_subject_rffrac_dict[subject]
	if 'M' in subdict and 'B' in subdict:
		backhed_rffrac_diff[0].append(subdict['B']-subdict['M'])
	if 'B' in subdict and '4M' in subdict:
		backhed_rffrac_diff[1].append(subdict['4M']-subdict['B'])
	if '4M' in subdict and '12M' in subdict:
		backhed_rffrac_diff[2].append(subdict['12M']-subdict['4M'])

# Ferretti
ferretti_subject_rffrac_dict = {}
for i in range(len(ferretti_ref_freq_frac)):
	short = ferretti_cohort_names[i]
	for sample in ferretti_ref_freq_frac[i].keys():
		subject = sample_subject_map[sample]
		if subject not in ferretti_subject_rffrac_dict:
			ferretti_subject_rffrac_dict[subject] = {}
		ferretti_subject_rffrac_dict[subject][short] = ferretti_ref_freq_frac[i][sample]

# 0: t0-t1, 1: t1-t2, 2: t2-t3, 3: t3-t4, 4: t4-t5
ferretti_rffrac_diff = [[] for i in range(5)]

for subject in ferretti_subject_rffrac_dict:
	subdict = ferretti_subject_rffrac_dict[subject]
	if 'M0' in subdict and 'I1' in subdict:
		ferretti_rffrac_diff[0].append(subdict['I1']-subdict['M0'])
	if 'I1' in subdict and 'I2' in subdict:
		ferretti_rffrac_diff[1].append(subdict['I2']-subdict['I1'])
	if 'I2' in subdict and 'I3' in subdict:
		ferretti_rffrac_diff[2].append(subdict['I3']-subdict['I2'])
	if 'I3' in subdict and 'I4' in subdict:
		ferretti_rffrac_diff[3].append(subdict['I4']-subdict['I3'])
	if 'I4' in subdict and 'I5' in subdict:
		ferretti_rffrac_diff[4].append(subdict['I5']-subdict['I4'])

# Yassour
yassour_subject_rffrac_dict = {}
for i in range(len(yassour_ref_freq_frac)):
	short = yassour_cohort_names[i]
	for sample in yassour_ref_freq_frac[i].keys():
		subject = sample_subject_map[sample]
		if subject not in yassour_subject_rffrac_dict:
			yassour_subject_rffrac_dict[subject] = {}
		yassour_subject_rffrac_dict[subject][short] = yassour_ref_freq_frac[i][sample]

# Yassour: M:delivery-C:birth, C:birth-C:week2, C:week2-C:month1, C:month1-C:month2, C:month1-C:month3
yassour_rffrac_diff = [[] for i in range(5)]

for subject in yassour_subject_rffrac_dict:
	subdict = yassour_subject_rffrac_dict[subject]
	if 'MBirth' in subdict and 'CBirth' in subdict:
		yassour_rffrac_diff[0].append(subdict['CBirth']-subdict['MBirth'])
	if 'CBirth' in subdict and 'C14' in subdict:
		yassour_rffrac_diff[1].append(subdict['C14']-subdict['CBirth'])
	if 'C14' in subdict and 'C1' in subdict:
		yassour_rffrac_diff[2].append(subdict['C1']-subdict['C14'])
	if 'C1' in subdict and 'C2' in subdict:
		yassour_rffrac_diff[3].append(subdict['C2']-subdict['C1'])
	if 'C2' in subdict and 'C3' in subdict:
		yassour_rffrac_diff[4].append(subdict['C3']-subdict['C2'])

###################################
#
# Boxplots for within-subject differences between cohorts
# Backhed: mother-infant, birth-4month, 4month-12month
# Ferretti: t0-t1, t1-t2, t2-t3, t3-t4, t4-t5
# Yassour: M:delivery-C:birth, C:birth-C:week2, C:week2-C:month1, C:month1-C:month2, C:month2-C:month3
#
###################################

plot_data = backhed_rffrac_diff
plot_names = ["M-birth", "birth-m4", "m4-m12"]
plot_names_n = ["{}\nn={}".format(a_, b_) for a_, b_ in zip(plot_names, [len(x) for x in plot_data])]

ax4.axhline(y=0, color='#e3e3e3', linestyle='--')
ax4.boxplot(plot_data)
ax4.set_ylabel("Within-subject\nchange in polymorphism")
ax4.set_xticklabels(plot_names_n)

plot_data = ferretti_rffrac_diff
plot_names = ["M-d1","d1-d3", "d3-w1", "w1-m1", "m1-m4"]
plot_names_n = ["{}\nn={}".format(a_, b_) for a_, b_ in zip(plot_names, [len(x) for x in plot_data])]

ax5.axhline(y=0, color='#e3e3e3', linestyle='--')
ax5.boxplot(plot_data)
ax5.set_xticklabels(plot_names_n)

plot_data = yassour_rffrac_diff
plot_names = ["M-birth", "birth-w2", "w2-m1","m1-m2","m2-m3"]
plot_names_n = ["{}\nn={}".format(a_, b_) for a_, b_ in zip(plot_names, [len(x) for x in plot_data])]

ax6.axhline(y=0, color='#e3e3e3', linestyle='--')
ax6.boxplot(plot_data)
ax6.set_xticklabels(plot_names_n)

####
#
# Save figure
#
####

plt.tight_layout()
sys.stderr.write("Saving figure...\t")
fig.savefig('%s/%s_%s_polymorphism_cohort_comparison.pdf' % (parse_midas_data.analysis_directory, upper_threshold, species_name), bbox_inches='tight')
sys.stderr.write("Done!\n")
