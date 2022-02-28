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
import os.path


minfreq=0.8
upper_threshold = minfreq
lower_threshold = 1 - upper_threshold
min_coverage = 20

good_species_list=parse_midas_data.parse_good_species_list()

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sys.stderr.write("Done!\n") 


def set_box_color(bp, color):
		plt.setp(bp['boxes'], color=color)
		plt.setp(bp['whiskers'], color=color, linestyle='-')
		plt.setp(bp['caps'], color=color)
		plt.setp(bp['medians'], color=color)		
		plt.setp(bp['fliers'], color=color, markersize=3)


# set up dataframe to store data
dataset_cohort={}
dataset_cohort['backhed']=['B','4M','12M','M']
dataset_cohort['yassour']=['MGest', 'MBirth', 'M3', 'CBirth', 'C14', 'C1', 'C2', 'C3']
dataset_cohort['ferretti']=['M0','I1','I2','I3','I4','I5']

ref_freq_frac_QP={}
for dataset in dataset_cohort:		 
		for age in dataset_cohort[dataset]:		
				ref_freq_frac_QP[age]=numpy.array([])

ref_freq_frac_nonQP={}
for dataset in dataset_cohort:		 
		for age in dataset_cohort[dataset]:		
				ref_freq_frac_nonQP[age]=numpy.array([])

for species_name in good_species_list:
	
	# Load SNP information for species_name
	sys.stderr.write("Loading SFSs for %s...\t" % species_name)
	_, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
	sys.stderr.write("Done!\n")
	
	# Load genomic coverage distributions
	sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
	median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
	sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
	samples = numpy.array(samples)
	
	# Get QP and non-QP samples lists
	samples_QP = diversity_utils.calculate_haploid_samples(species_name)
	samples_nonQP = numpy.setdiff1d(samples, samples_QP)
	
	# Only plot samples above a certain depth threshold
	desired_samples = samples[(median_coverages>=min_coverage)]
	desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])
	
	if len(desired_samples) <= 0:
		continue
	
	# Calculate within polymorphism rates	
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
		#within_sites, between_sites, total_sites = sfs_utils.calculate_singleton_rates_from_sfs_map(sfs_map[sample], lower_threshold=0, upper_threshold=0.1)
		if total_sites <= 0:
			continue
		
		within_rate = within_sites*1.0/total_sites
		between_rate = between_sites*1.0/total_sites
		between_rates.append(between_rate)
		within_rates.append(within_rate)
		#
		# Calculate 95% confidence intervals
		within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05)
		within_rate_lowers.append(within_rate_lower)
		within_rate_uppers.append(within_rate_upper)
		# 
		depths, counts = sfs_utils.calculate_depth_distribution_from_sfs_map(sfs_map[sample])
		dlower, dupper = stats_utils.calculate_IQR_from_distribution(depths, counts)
		dmedian = stats_utils.calculate_median_from_distribution(depths,counts)
		#
		depth_lowers.append(dlower)
		depth_uppers.append(dupper)
		median_depths.append(dmedian)
		sample_names.append(sample)
	
	# Sort them all in descending order of within-host diversity		
	within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths, depth_lowers, depth_uppers, sample_names = (numpy.array(x) for x in zip(*sorted(zip(within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths,depth_lowers, depth_uppers, sample_names),reverse=True)))
	
	within_rate_lowers = numpy.clip(within_rate_lowers, 1e-09,1)
	
	#store refreq for each age for QP samples:
	for dataset in dataset_cohort:
		for age in dataset_cohort[dataset]:
			sample_names_cohort = sample_utils.get_sample_names(dataset, age)
			for sample_name in sample_names_cohort:
				idx=numpy.where(sample_names==sample_name)
				if sample_name in samples_QP:
					ref_freq_frac_QP[age] = numpy.append(ref_freq_frac_QP[age], within_rates[idx])
				elif sample_name in samples_nonQP:
					ref_freq_frac_nonQP[age] = numpy.append(ref_freq_frac_nonQP[age], within_rates[idx])

ref_freq_frac_QP_pickle_fn = "%s/pickles/ref_freq_frac_QP.pkl" % config.data_directory
ref_freq_frac_nonQP_pickle_fn = "%s/pickles/ref_freq_frac_nonQP.pkl" % config.data_directory

import pickle

pickle.dump(ref_freq_frac_QP, open(ref_freq_frac_QP_pickle_fn, 'wb'))
pickle.dump(ref_freq_frac_nonQP, open(ref_freq_frac_nonQP_pickle_fn, 'wb'))

def plot_legend():

		plt.plot([0,.1],[0.9,0.9], color='#fed976',linewidth=2.0)
		plt.text(0.2, 0.9, 'Birth', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.85,0.85], color='#feb24c',linewidth=2.0) 
		plt.text(0.2, 0.85, '1 day', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.8,0.8], color='#fd8d3c',linewidth=2.0) 
		plt.text(0.2, 0.8, '3 days', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.75,0.75], color='#fc4e2a',linewidth=2.0) 
		plt.text(0.2, 0.75, '4 days', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.7,0.7], color='#e31a1c',linewidth=2.0) 
		plt.text(0.2, 0.7, '7 days', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.65,0.65], color='#b10026',linewidth=2.0) 
		plt.text(0.2, 0.65, '14 days', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.6,0.6], color='#bdd7e7',linewidth=2.0) 
		plt.text(0.2, 0.6, '1 month', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.55,0.55], color='#6baed6',linewidth=2.0) 
		plt.text(0.2, 0.55, '2 months', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.5,0.5], color='#08519c',linewidth=2.0) 
		plt.text(0.2, 0.5, '3 months', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.45,0.45], color='#08306b',linewidth=2.0) 
		plt.text(0.2, 0.45, '4 months', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.4,0.4], color='#addd8e',linewidth=2.0) 
		plt.text(0.2, 0.4, '12 months', horizontalalignment='center',verticalalignment='center')
		plt.plot([0,.1],[0.35,0.35], color='pink',linewidth=2.0) 
		plt.text(0.2, 0.35, 'mothers', horizontalalignment='center',verticalalignment='center')
		
		plt.xlim(0,1)


######################################################
# after iterating through all species, now plot
######################################################

# combine all 1 month samples
ref_freq_frac_QP['1_month']=numpy.concatenate((ref_freq_frac_QP['C1'],ref_freq_frac_QP['I4']), axis=0)
ref_freq_frac_nonQP['1_month']=numpy.concatenate((ref_freq_frac_nonQP['C1'],ref_freq_frac_nonQP['I4']), axis=0)

# combine all 4 month sample
ref_freq_frac_QP['4_months']=numpy.concatenate((ref_freq_frac_QP['4M'],ref_freq_frac_QP['I5']), axis=0)
ref_freq_frac_nonQP['4_months']=numpy.concatenate((ref_freq_frac_nonQP['4M'],ref_freq_frac_nonQP['I5']), axis=0)

# combine all mother samples
ref_freq_frac_QP['all_mothers']=numpy.concatenate((ref_freq_frac_QP['M'],ref_freq_frac_QP['MBirth'],ref_freq_frac_QP['M0']), axis=0)
ref_freq_frac_nonQP['all_mothers']=numpy.concatenate((ref_freq_frac_nonQP['M'],ref_freq_frac_nonQP['MBirth'],ref_freq_frac_nonQP['M0']), axis=0)

age_order=['CBirth','I1','I2','B','I3', 'C14', '1_month','C2','C3','4_months','12M','all_mothers']


colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#08306b', '#addd8e','pink']

###
# Boxplot of just QP samples

matplotlib.rcParams.update({'font.size': 9})
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 4

plt.subplot(1,2, 1)

col_i=0
for age in age_order:
					print age + '\t' + str(len(ref_freq_frac_QP[age])) + '\t' + str(len(ref_freq_frac_nonQP[age])) + '\t' + str(len(ref_freq_frac_QP[age])*1.0/(len(ref_freq_frac_QP[age]) + len(ref_freq_frac_nonQP[age]))) 
					bp = plt.boxplot(ref_freq_frac_QP[age], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
					set_box_color(bp, colors[col_i]) 
					col_i+=1

# dummy boxplot to get the colors to work for last one
bp = plt.boxplot(ref_freq_frac_QP[age], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
#patch.set_facecolor('white')
set_box_color(bp, 'white') 
plt.xlim(0,len(age_order)+2)
plt.yscale('log')

ax=plt.subplot(1, 2, 2)

plot_legend()
plt.xlim(0,1)
plt.tight_layout()
plt.axis('off')
plt.savefig(os.path.expanduser('~/mother_infant/analysis/within_host_piS_QP.png'), dpi=600)
plt.close()				 

###
# Boxplot of just nonQP samples

matplotlib.rcParams.update({'font.size': 9})
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 4

plt.subplot(1,2, 1)

col_i=0
for age in age_order:
					bp = plt.boxplot(ref_freq_frac_nonQP[age], positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
					set_box_color(bp, colors[col_i]) 
					col_i+=1

# dummy boxplot to get the colors to work for last one
bp = plt.boxplot(ref_freq_frac_QP[age], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
#patch.set_facecolor('white')
set_box_color(bp, 'white') 
plt.xlim(0,len(age_order)+2)
plt.yscale('log')

ax=plt.subplot(1, 2, 2)

plot_legend()
plt.xlim(0,1)
plt.tight_layout()
plt.axis('off')
plt.tight_layout()
plt.savefig(os.path.expanduser('~/mother_infant/analysis/within_host_piS_nonQP.png'), dpi=600)
plt.close()				 



###
# Boxplot of all samples

matplotlib.rcParams.update({'font.size': 9})
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 8
fig_size[1] = 4

plt.subplot(1,2, 1)

col_i=0
for age in age_order:
					bp = plt.boxplot(numpy.append(ref_freq_frac_nonQP[age], ref_freq_frac_QP[age]), positions=[col_i+2], sym='.',patch_artist=True, widths=0.5, labels=[''])
					set_box_color(bp, colors[col_i]) 
					col_i+=1

# dummy boxplot to get the colors to work for last one
bp = plt.boxplot(ref_freq_frac_QP[age], positions=[col_i+10], sym='.',patch_artist=True, widths=0.5, labels=[''])
set_box_color(bp, 'white') 
plt.xlim(0,len(age_order)+2)
plt.yscale('log')


ax=plt.subplot(1, 2, 2)

plot_legend()
plt.xlim(0,1)
plt.tight_layout()
plt.axis('off')
plt.savefig(os.path.expanduser('~/mother_infant/analysis/within_host_piS_all.png'), dpi=600)
plt.close()				 


###################
# Try histogram plot #

for age in	age_order: 
		bin_defs=numpy.arange(0,0.02, step=0.0001)
		n, bins, patches = plt.hist(x=ref_freq_frac_QP[age], color='#0504aa', bins=bin_defs, alpha=0.5)
		n, bins, patches = plt.hist(x=ref_freq_frac_nonQP[age], color='#fed976', bins=bin_defs, alpha=0.5)


		plt.grid(axis='y', alpha=0.75)
		plt.xlabel('Value')
		plt.ylabel('Frequency')
		plt.title('Within host polymorphism levels %s' %age)
		plt.text(23, 45, r'$\mu=15, b=3$')
		maxfreq = n.max()
		plt.tight_layout()
		plt.savefig(os.path.expanduser('~/mother_infant/analysis/within_host_piS_hist_QP_%s.png' %age), dpi=600)
		plt.close()


sys.stderr.write("Done!\n")
