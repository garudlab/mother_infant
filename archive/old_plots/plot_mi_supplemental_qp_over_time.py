import matplotlib	 
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import sample_utils
import os.path
import pylab
import sys
import numpy as np
import diversity_utils
import gene_diversity_utils
import calculate_temporal_changes
import calculate_substitution_rates
import stats_utils
import sfs_utils
from collections import defaultdict

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, random, choice, multinomial
import matplotlib.colors as mcolors # Redundant?
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']	= False
mpl.rcParams['legend.fontsize']	 = 'small'

species_name = "Bacteroides_vulgatus_57955"

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--modification-threshold", type=int, help="max number of SNV differences before calling a modification", default=config.modification_difference_threshold)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize
modification_difference_threshold = args.modification_threshold
replacement_difference_threshold = config.replacement_difference_threshold
twin_modification_difference_threshold = config.twin_modification_difference_threshold
twin_replacement_difference_threshold = config.twin_replacement_difference_threshold

twin_modification_difference_threshold = 1e06
twin_replacement_difference_threshold = 1e06

################################################################################

#####################
#
# Settings for calculation:
#
#####################

min_coverage = config.min_median_coverage
min_sample_size = 3
min_haploid_sample_size = 10

#####
#
# Settings for different cohorts we are looking at 
#
#####
cohorts = ["backhed", "ferretti", "yassour"]
countries = ["Sweden", "Italy", "Finland"]
country_cohort_map = {country: cohort for country,cohort in zip(countries,cohorts)}

# Change later
modification_difference_thresholds = {"backhed": modification_difference_threshold, "ferretti": modification_difference_threshold, "yassour": modification_difference_threshold}

replacement_difference_thresholds = {"backhed": replacement_difference_threshold, "ferretti": replacement_difference_threshold, "yassour": replacement_difference_threshold}

################################
#
# Set up figures
#
################################


####################################################
#
# Set up Supplemental Fig (temporal haploid classification)
#
####################################################
# This figure spreads them all out

pylab.figure(6,figsize=(7.5,6))
fig6 = pylab.gcf()
# make three panels panels
outer_grid6	 = gridspec.GridSpec(1,3,width_ratios=[1,1,1],wspace=0.2)

backhed_haploid_axis = plt.Subplot(fig6, outer_grid6[0])
fig6.add_subplot(backhed_haploid_axis)
backhed_haploid_axis.set_xlabel('Backhed timepoint pairs')

ferretti_haploid_axis = plt.Subplot(fig6, outer_grid6[1])
fig6.add_subplot(ferretti_haploid_axis)
ferretti_haploid_axis.set_xlabel('Ferretti timepoint pairs')

yassour_haploid_axis = plt.Subplot(fig6, outer_grid6[2])
fig6.add_subplot(yassour_haploid_axis)
yassour_haploid_axis.set_xlabel('Yassour timepoint pairs')

################################
#
# Now do calculation
#
################################

backhed_species_qp_counts = {}
ferretti_species_qp_counts = {}
yassour_species_qp_counts = {}

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sample_country_map = sample_utils.parse_sample_country_map()
sys.stderr.write("Done!\n")

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
		good_species_list = ["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318"]

num_passed_species = 0

all_samples = sample_order_map.keys()
same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, all_samples)

for species_name in good_species_list:
		sys.stderr.write("\nProcessing %s...\n" % species_name)
		
		# First we have to enumerate QP pairs in each cohort
		sys.stderr.write("Enumerating QP pairs...\n")
		
		# list of samples that meet coverage criteria for this species
		highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
		
		# list of samples that meet QP criteria for this species
		haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))
		
		#print len(all_samples), len(highcoverage_samples), len(haploid_samples)
		
		if len(haploid_samples) < min_haploid_sample_size:
				continue
		
		#all_samples = list(haploid_samples)
		
		hmp_sample_size = 0
		
		qp_sample_sets = {cohort: set() for cohort in cohorts}
		qp_counts = {cohort:[0,0,0,0] for cohort in cohorts}
		
		for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):		
				i = same_subject_idxs[0][sample_pair_idx]
				j = same_subject_idxs[1][sample_pair_idx]
				
				sample_i = all_samples[i]
				sample_j = all_samples[j]
				
				country = sample_country_map[sample_i]
				
				if country not in countries:
						continue
				
				# Figure out cohort
				cohort = country_cohort_map[country]
				
				# Figure out QP status of pair				
				if not ((sample_i in highcoverage_samples) and (sample_j in highcoverage_samples)):
						# NOT both are highcoverage samples						
						if ((sample_i in highcoverage_samples) or (sample_j in highcoverage_samples)):
								# One sample is high coverage
								qp_counts[cohort][0] += 1
						else:
								# Neither sample is high coverage, ignore
								pass						
				else:						
						# Both are highcoverage samples						
						if (sample_i in haploid_samples) and (sample_j in haploid_samples):								
								# Both are QP samples!								
								qp_counts[cohort][1] += 1
								qp_sample_sets[cohort].add(sample_i)
								qp_sample_sets[cohort].add(sample_j) 
								#print sample_i, sample_j						
						elif (sample_i not in haploid_samples) and (sample_j not in haploid_samples):
								# pair that is non-QP at both timepoints
								qp_counts[cohort][2] += 1						
						else:
								# pair that is QP at one timepoint and non-QP at another
								qp_counts[cohort][3] += 1		
		sys.stderr.write("Done!\n")
		
		for cohort in cohorts:		
				print ("%s:" % cohort), qp_counts[cohort][1], "QP pairs,", qp_counts[cohort][2], "non-QP pairs,", qp_counts[cohort][3], "mixed pairs", qp_counts[cohort][0], "species dropouts."
		
		combined_sample_set = set()
		for cohort in cohorts:
				combined_sample_set.update(qp_sample_sets[cohort])
		combined_samples = list(sorted(combined_sample_set))
		combined_sample_idx_map = {combined_samples[i] : i for i in xrange(0,len(combined_samples))}
		
		qp_sample_lists = {cohort: list(sorted(qp_sample_sets[cohort])) for cohort in cohorts}
		
		backhed_sample_size = len(qp_sample_sets['backhed'])
		
		if backhed_sample_size < min_sample_size:
				continue
		
		backhed_species_qp_counts[species_name] = qp_counts['backhed']
		ferretti_species_qp_counts[species_name] = qp_counts['ferretti']
		yassour_species_qp_counts[species_name] = qp_counts['yassour']
		
		sys.stderr.write("Proceeding with %d HMP longitudinal comparisons!\n" % (backhed_sample_size))

### Now plot temporal qp figures

species_names = backhed_species_qp_counts.keys()
species_names = list(sorted(species_names, key=lambda s: sum(backhed_species_qp_counts[s])))

ys = np.arange(0,len(species_names))

yticklabels = []

for y,species_name in zip(ys,species_names):
		yticklabels.append(species_name)
		
		total_samples = sum(backhed_species_qp_counts[species_name])
		
		if total_samples>0:		
				qp_samples = backhed_species_qp_counts[species_name][1]
				non_qp_samples = backhed_species_qp_counts[species_name][2]
				mixed_samples = backhed_species_qp_counts[species_name][3]
				dropout_samples = backhed_species_qp_counts[species_name][0]
				
				backhed_haploid_axis.barh([y],[qp_samples],linewidth=0, color='#08519c')
				
				backhed_haploid_axis.barh([y],[non_qp_samples], left=[qp_samples], linewidth=0, color='#de2d26')
				
				backhed_haploid_axis.barh([y],[mixed_samples], left=[qp_samples+non_qp_samples], linewidth=0, color='#8856a7')
				
				backhed_haploid_axis.barh([y],[dropout_samples], left=[qp_samples+non_qp_samples+mixed_samples], linewidth=0, color='0.7')
		
		total_samples = sum(ferretti_species_qp_counts[species_name])
		
		if total_samples>0:		
				qp_samples = ferretti_species_qp_counts[species_name][1]
				non_qp_samples = ferretti_species_qp_counts[species_name][2]
				mixed_samples = ferretti_species_qp_counts[species_name][3]
				dropout_samples = ferretti_species_qp_counts[species_name][0]
				
				ferretti_haploid_axis.barh([y],[qp_samples],linewidth=0,label='QP->QP', color='#08519c')
				
				ferretti_haploid_axis.barh([y],[non_qp_samples], left=[qp_samples], linewidth=0, color='#de2d26')
				
				ferretti_haploid_axis.barh([y],[mixed_samples], left=[qp_samples+non_qp_samples], linewidth=0, color='#8856a7')
				
				ferretti_haploid_axis.barh([y],[dropout_samples], left=[qp_samples+non_qp_samples+mixed_samples], linewidth=0, color='0.7')
		
		total_samples = sum(yassour_species_qp_counts[species_name])
		
		if total_samples>0:		
				qp_samples = yassour_species_qp_counts[species_name][1]
				non_qp_samples = yassour_species_qp_counts[species_name][2]
				mixed_samples = yassour_species_qp_counts[species_name][3]
				dropout_samples = yassour_species_qp_counts[species_name][0]
				
				yassour_haploid_axis.barh([y],[qp_samples],linewidth=0,label='QP->QP', color='#08519c')
				
				yassour_haploid_axis.barh([y],[non_qp_samples], left=[qp_samples], linewidth=0, color='#de2d26')
				
				yassour_haploid_axis.barh([y],[mixed_samples], left=[qp_samples+non_qp_samples], linewidth=0, color='#8856a7')
				
				yassour_haploid_axis.barh([y],[dropout_samples], left=[qp_samples+non_qp_samples+mixed_samples], linewidth=0, color='0.7')

backhed_haploid_axis.yaxis.tick_left()
backhed_haploid_axis.xaxis.tick_bottom()
	
ferretti_haploid_axis.yaxis.tick_left()
ferretti_haploid_axis.xaxis.tick_bottom()

yassour_haploid_axis.yaxis.tick_left()
yassour_haploid_axis.xaxis.tick_bottom()

backhed_haploid_axis.set_yticks(ys+0.5)
ferretti_haploid_axis.set_yticks(ys+0.5)
yassour_haploid_axis.set_yticks(ys+0.5)

backhed_haploid_axis.set_ylim([-1,len(ys)])
ferretti_haploid_axis.set_ylim([-1,len(ys)])
yassour_haploid_axis.set_ylim([-1,len(ys)])

backhed_haploid_axis.set_yticklabels(yticklabels,fontsize=5)
ferretti_haploid_axis.set_yticklabels([])
yassour_haploid_axis.set_yticklabels([])

backhed_haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)
ferretti_haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)
yassour_haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)

backhed_haploid_axis.set_xlim([0,250])
ferretti_haploid_axis.set_xlim([0,50])
yassour_haploid_axis.set_xlim([0,150])

### Do stuff for legend
backhed_haploid_axis.barh([-10],[1],linewidth=0,label='QP->QP', color='#08519c')
backhed_haploid_axis.barh([-10],[1],linewidth=0,label='non->non', color='#de2d26')
backhed_haploid_axis.barh([-10],[1],linewidth=0,label='mixed', color='#8856a7')
backhed_haploid_axis.barh([-10],[1],linewidth=0,label='dropout', color='0.7')
backhed_haploid_axis.legend(loc='lower right',frameon=False)

###################################
#
# Set up next figure (23 panels, 23 rows / 1 col)
#
###################################

# Use the same 23 most prevalent species from above
# species_names = ['Bacteroides_fragilis_56548', 'Bacteroides_finegoldii_57739', 'Sutterella_wadsworthensis_56828', 'Bifidobacterium_catenulatum_58257', 'Bacteroides_xylanisolvens_57185', 'Bacteroides_faecis_58503', 'Bacteroides_thetaiotaomicron_56941', 'Bacteroides_stercoris_56735', 'Parabacteroides_merdae_56972', 'Bifidobacterium_pseudocatenulatum_57754', 'Bacteroides_caccae_53434', 'Enterococcus_faecalis_56297', 'Bifidobacterium_breve_57133', 'Bifidobacterium_bifidum_55065', 'Bifidobacterium_adolescentis_56815', 'Bacteroides_ovatus_58035', 'Ruminococcus_gnavus_57638', 'Bacteroides_fragilis_54507', 'Parabacteroides_distasonis_56985', 'Bacteroides_uniformis_57318', 'Bifidobacterium_longum_57796', 'Escherichia_coli_58110', 'Bacteroides_vulgatus_57955']
num_species = len(species_names)

figBT, axesBT = plt.subplots(num_species, 1, figsize=(10,12))
figFT, axesFT = plt.subplots(num_species, 1, figsize=(11,16))
figYT, axesYT = plt.subplots(num_species, 1, figsize=(11,16))

# Colormap
qp_colors = ['black', 'gray', 'orange', 'blue']

def gen_colormap(matrix):
	qp_colors = {-1: 'black', 0: 'gray', 1: 'orange', 2: 'blue'}
	use_colors = []
	for row in matrix:
		for val in [-1,0,1,2]:		
			if val in row:
				use_colors.append(qp_colors[val])
	cm = colors.LinearSegmentedColormap.from_list('qp_color_list', use_colors, N = len(use_colors))
	return cm

################################
#
# Individual subjects
#
################################

qp_cats = ['nonQP','QP','lowcov']

backhed_samples = sample_utils.get_sample_names('backhed', 'all')
backhed_cohorts = ['M','B','4M','12M']
nb = len(backhed_cohorts)

ferretti_samples = sample_utils.get_sample_names('ferretti', 'all')
ferretti_cohorts = ['M0','I1','I2','I3','I4','I5']
nf = len(ferretti_cohorts)

yassour_samples = sample_utils.get_sample_names('yassour', 'all')
yassour_cohorts = ['MGest', 'MBirth', 'M3', 'CBirth', 'C14', 'C1', 'C2', 'C3']
ny = len(yassour_cohorts)

btp = ['M1','I1','I2','I3']
backhed_subject_qp_agg = {i: {cat: defaultdict(int) for cat in qp_cats} for i in btp}
ftp = ['M1','I1','I2','I3','I4','I5']
ferretti_subject_qp_agg = {i: {cat: defaultdict(int) for cat in qp_cats} for i in ftp}
ytp = ['M1','M2','M3','C1','C2','C3','C4','C5']
yassour_subject_qp_agg = {i: {cat: defaultdict(int) for cat in qp_cats} for i in ytp}

for i in range(num_species):
	species_name = species_names[i]
	
	# list of samples that meet coverage criteria for this species
	highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
	
	# list of samples that meet QP criteria for this species
	haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))
	
	# Backhed	
	backhed_subject_qp = {}
	sys.stderr.write("\nProcessing Backhed %s...\n" % species_name)
	
	for sample in backhed_samples:
		subject, order = sample_order_map[sample]
		timept = subject[-1] + str(order)
		if subject not in backhed_subject_qp:
			backhed_subject_qp[subject] = { tp:-1 for tp in btp }		
		if sample in highcoverage_samples:
			if sample in haploid_samples: # QP
				backhed_subject_qp[subject][timept] = 1
				backhed_subject_qp_agg[timept]['QP'][species_name] += 1
			else: # Non-QP but high coverage
				backhed_subject_qp[subject][timept] = 2
				backhed_subject_qp_agg[timept]['nonQP'][species_name] += 1
		else: # Low coverage for this species
			backhed_subject_qp[subject][timept] = 0
			backhed_subject_qp_agg[timept]['lowcov'][species_name] += 1
	
	col_labels = backhed_subject_qp.keys() # Subject IDs
	row_labels = backhed_cohorts # Cohort names
	backhed_matrix = np.array([backhed_subject_qp[j].values() for j in col_labels])
	backhed_matrix.sort(axis=0)
	cm = gen_colormap(backhed_matrix)
	backhed_matrix = backhed_matrix.transpose()
	
	axesBT[i].matshow(backhed_matrix, cmap=cm)
	axesBT[i].set_xticks([], [])
	axesBT[i].set_yticklabels([''] + backhed_cohorts)
	axesBT[i].tick_params(axis='both', which='both',length=0)
	axesBT[i].yaxis.set_label_position('right')
	axesBT[i].grid(which='both', axis='x', linestyle='-', color='k', linewidth=0.5)
	axesBT[i].set_ylabel(species_name, rotation=0)
	axesBT[i].yaxis.set_label_coords(-0.17, 0.5)
	
	# Ferretti
	ferretti_subject_qp = {}
	sys.stderr.write("\nProcessing Ferretti %s...\n" % species_name)
	
	for sample in ferretti_samples:
		subject, order = sample_order_map[sample]
		timept = subject[-2] + str(order)
		if subject not in ferretti_subject_qp:
			ferretti_subject_qp[subject] = { tp:-1 for tp in ftp }		
		if sample in highcoverage_samples:
			if sample in haploid_samples: # QP
				ferretti_subject_qp[subject][timept] = 1
				ferretti_subject_qp_agg[timept]['QP'][species_name] += 1
			else: # Non-QP but high coverage
				ferretti_subject_qp[subject][timept] = 2
				ferretti_subject_qp_agg[timept]['nonQP'][species_name] += 1
		else: # Low coverage for this species
			ferretti_subject_qp[subject][timept] = 0
			ferretti_subject_qp_agg[timept]['lowcov'][species_name] += 1
	
	col_labels = ferretti_subject_qp.keys() # Subject IDs
	row_labels = ferretti_cohorts # Cohort names
	ferretti_matrix = np.array([ferretti_subject_qp[j].values() for j in col_labels])
	ferretti_matrix.sort(axis=0)
	cm = gen_colormap(ferretti_matrix)
	ferretti_matrix = ferretti_matrix.transpose()
	
	axesFT[i].matshow(ferretti_matrix, cmap=cm)
	axesFT[i].set_xticks([], [])
	axesFT[i].set_yticklabels([''] + ferretti_cohorts)
	axesFT[i].tick_params(axis='both', which='both',length=0)
	axesFT[i].yaxis.set_label_position('right')
	axesFT[i].grid(which='both', axis='x', linestyle='-', color='k', linewidth=0.5)
	axesFT[i].set_ylabel(species_name, rotation=0)
	axesFT[i].yaxis.set_label_coords(-0.22, 0.5)
	
	# Yassour
	yassour_subject_qp = {}
	sys.stderr.write("\nProcessing Yassour %s...\n" % species_name)
	
	for sample in yassour_samples:
		subject, order = sample_order_map[sample]
		timept = subject[-1] + str(order)
		if subject not in yassour_subject_qp:
			yassour_subject_qp[subject] = { tp:-1 for tp in ytp }		
		if sample in highcoverage_samples:
			if sample in haploid_samples: # QP
				yassour_subject_qp[subject][timept] = 1
				yassour_subject_qp_agg[timept]['QP'][species_name] += 1
			else: # Non-QP but high coverage
				yassour_subject_qp[subject][timept] = 2
				yassour_subject_qp_agg[timept]['nonQP'][species_name] += 1
		else: # Low coverage for this species
			yassour_subject_qp[subject][timept] = 0
			yassour_subject_qp_agg[timept]['lowcov'][species_name] += 1
	
	col_labels = yassour_subject_qp.keys() # Subject IDs
	row_labels = yassour_cohorts # Cohort names
	yassour_matrix = np.array([yassour_subject_qp[j].values() for j in col_labels])
	yassour_matrix.sort(axis=0)
	cm = gen_colormap(yassour_matrix)
	yassour_matrix = yassour_matrix.transpose()
	
	axesYT[i].matshow(yassour_matrix, cmap=cm)
	axesYT[i].set_xticks([], [])
	axesYT[i].set_yticklabels([''] + yassour_cohorts)
	axesYT[i].tick_params(axis='both', which='both',length=0)
	axesYT[i].yaxis.set_label_position('right')
	axesYT[i].grid(which='both', axis='x', linestyle='-', color='k', linewidth=0.5)
	axesYT[i].set_ylabel(species_name, rotation=0)
	axesYT[i].yaxis.set_label_coords(-0.25, 0.5)

legend_labels = ["No sample","Low coverage","QP","Non-QP"]
legend_patches = [mpatches.Patch(color=qp_color, label=[label]) for qp_color, label in zip(qp_colors, legend_labels)]

figBT.suptitle("QP Changes for Backhed samples")
axesBT[0].set_xlabel("Subjects (n=98)")
axesBT[0].xaxis.set_label_position('top')
axesBT[0].legend(handles=legend_patches, bbox_to_anchor=(1.13,1),loc="upper right")

figFT.suptitle("QP Changes for Ferretti samples")
axesFT[0].set_xlabel("Subjects")
axesFT[0].xaxis.set_label_position('top')
axesFT[0].legend(handles=legend_patches, bbox_to_anchor=(1.15,1),loc="upper right")

figYT.suptitle("QP Changes for Yassour samples")
axesYT[0].set_xlabel("Subjects")
axesYT[0].xaxis.set_label_position('top')
axesYT[0].legend(handles=legend_patches, bbox_to_anchor=(1.15,1),loc="upper right")

###################################
#
# Set up next figure: Number of QP hosts by cohort
#
###################################

figAgg1, axesAgg1 = plt.subplots(1, 3, sharey='row', figsize=(12,4),gridspec_kw={'width_ratios': [2,3,4], 'wspace': 0.02})

###################################
#
# Aggregate over all species
# Separate by dataset
# Just show QP vs non-QP (both high coverage) for each cohort
#
###################################

agg_count_QP = {'backhed': {i:0 for i in btp} ,'ferretti': {i:0 for i in ftp}, 'yassour': {i:0 for i in ytp}}
agg_count_nonQP = {'backhed': {i:0 for i in btp} ,'ferretti': {i:0 for i in ftp}, 'yassour': {i:0 for i in ytp}}
agg_count_lowcov = {'backhed': {i:0 for i in btp} ,'ferretti': {i:0 for i in ftp}, 'yassour': {i:0 for i in ytp}}

# Backhed
for tp in backhed_subject_qp_agg.keys():
	for count in backhed_subject_qp_agg[tp]['QP'].values():
		agg_count_QP['backhed'][tp] += count
	for count in backhed_subject_qp_agg[tp]['nonQP'].values():
		agg_count_nonQP['backhed'][tp] += count
	for count in backhed_subject_qp_agg[tp]['lowcov'].values():
		agg_count_lowcov['backhed'][tp] += count

# Ferretti
for tp in ferretti_subject_qp_agg.keys():
	for count in ferretti_subject_qp_agg[tp]['QP'].values():
		agg_count_QP['ferretti'][tp] += count
	for count in ferretti_subject_qp_agg[tp]['nonQP'].values():
		agg_count_nonQP['ferretti'][tp] += count
	for count in ferretti_subject_qp_agg[tp]['lowcov'].values():
		agg_count_lowcov['ferretti'][tp] += count

# Yassour
for tp in yassour_subject_qp_agg.keys():
	for count in yassour_subject_qp_agg[tp]['QP'].values():
		agg_count_QP['yassour'][tp] += count
	for count in yassour_subject_qp_agg[tp]['nonQP'].values():
		agg_count_nonQP['yassour'][tp] += count
	for count in yassour_subject_qp_agg[tp]['lowcov'].values():
		agg_count_lowcov['yassour'][tp] += count

###################################
#
# Plot number QP samples (hosts) by timepoint cohort
#
###################################

barWidth = 0.2

# Backhed
barPos_lowcov = np.arange(nb) + barWidth
barPos_QP = [x + barWidth for x in barPos_lowcov]
barPos_nonQP = [x + barWidth for x in barPos_QP]
tps = ['M1','I1','I2','I3']
backhed_cohorts = ["Mothers\nDelivery", "3-5 Days","4 Months", "12 Months"]

bars = [agg_count_lowcov['backhed'][tp] for tp in tps]
axesAgg1[0].bar(barPos_lowcov, bars, color='gray', width=barWidth, edgecolor='white', label='Low coverage')

bars = [agg_count_QP['backhed'][tp] for tp in tps]
axesAgg1[0].bar(barPos_QP, bars, color='orange', width=barWidth, edgecolor='white', label='QP')

bars = [agg_count_nonQP['backhed'][tp] for tp in tps]
axesAgg1[0].bar(barPos_nonQP, bars, color='blue', width=barWidth, edgecolor='white', label='Non-QP')

axesAgg1[0].set_title("Backhed", fontweight='bold')
axesAgg1[0].set_ylabel("Number of host-species pairs\n(pooled across 23 prevalent species)")
axesAgg1[0].set_yscale('log')
axesAgg1[0].set_xticks(barPos_QP)
axesAgg1[0].set_xticklabels(backhed_cohorts)
axesAgg1[0].legend()

# Ferretti
barPos_lowcov = np.arange(nf) + barWidth
barPos_QP = [x + barWidth for x in barPos_lowcov]
barPos_nonQP = [x + barWidth for x in barPos_QP]
tps = ['M1','I1','I2','I3','I4','I5']
ferretti_cohorts = ["Mothers\nDelivery","1 Day", "3 Days", "1 Week", "1 Month", "4 Months"]

bars = [agg_count_lowcov['ferretti'][tp] for tp in tps]
axesAgg1[1].bar(barPos_lowcov, bars, color='gray', width=barWidth, edgecolor='white', label='Low coverage')

bars = [agg_count_QP['ferretti'][tp] for tp in tps]
axesAgg1[1].bar(barPos_QP, bars, color='orange', width=barWidth, edgecolor='white', label='QP')

bars = [agg_count_nonQP['ferretti'][tp] for tp in tps]
axesAgg1[1].bar(barPos_nonQP, bars, color='blue', width=barWidth, edgecolor='white', label='Non-QP')

axesAgg1[1].set_title("Ferretti", fontweight='bold')
# axesAgg1[1].set_ylabel("Number of hosts")
axesAgg1[1].set_yscale('log')
axesAgg1[1].set_xticks(barPos_QP)
axesAgg1[1].set_xticklabels(ferretti_cohorts)
axesAgg1[1].legend()

# Yassour
barPos_lowcov = np.arange(ny) + barWidth
barPos_QP = [x + barWidth for x in barPos_lowcov]
barPos_nonQP = [x + barWidth for x in barPos_QP]
tps = ['M1','M2','M3','C1','C2','C3','C4','C5']
yassour_cohorts = ["Mothers\nGestation","Mothers\nDelivery", "Mothers\n3 Months", "Birth", "2 Weeks", "1 Month", "2 Months", "3 Months"]

bars = [agg_count_lowcov['yassour'][tp] for tp in tps]
axesAgg1[2].bar(barPos_lowcov, bars, color='gray', width=barWidth, edgecolor='white', label='Low coverage')

bars = [agg_count_QP['yassour'][tp] for tp in tps]
axesAgg1[2].bar(barPos_QP, bars, color='orange', width=barWidth, edgecolor='white', label='QP')

bars = [agg_count_nonQP['yassour'][tp] for tp in tps]
axesAgg1[2].bar(barPos_nonQP, bars, color='blue', width=barWidth, edgecolor='white', label='Non-QP')

axesAgg1[2].set_title("Yassour", fontweight='bold')
# axesAgg1[2].set_ylabel("Number of hosts")
axesAgg1[2].set_yscale('log')
axesAgg1[2].set_xticks(barPos_QP)
axesAgg1[2].set_xticklabels(yassour_cohorts)
axesAgg1[2].legend()

###################################
#
# Save figures
#
###################################

sys.stderr.write("Saving figures...\t")
fig6.savefig('%s/temporal_haploid_mi.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

# figBT.savefig('%s/temporal_QP_by_subject_mi_backhed.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
# figFT.savefig('%s/temporal_QP_by_subject_mi_ferretti.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
# figYT.savefig('%s/temporal_QP_by_subject_mi_yassour.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

# figAgg1.savefig('%s/QP_by_cohort_agg_species_mi.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
