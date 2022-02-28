import matplotlib	 
matplotlib.use('Agg') 
import sample_utils as su
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
import diversity_utils
import gene_diversity_utils
import calculate_temporal_changes
import calculate_substitution_rates
import stats_utils
import sfs_utils
from collections import defaultdict
from numpy.random import choice

import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, random, choice, multinomial, shuffle
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.patches import Patch

import pickle
from os import path

# Plot utilities

# Label format function
def format_tp_label(cohort, tp_pair):
	tpa, tpb = tp_pair
	if tpa[0] == 'M' and tpb[0] == 'I' or tpa[0] == 'I' and tpb[0] == 'M':
		tp1, tp2 = (tpa, tpb) if tpa[0] == 'M' else (tpb, tpa)
	else:
		tp1, tp2 = (tpa, tpb) if (tpa < tpb) else (tpb, tpa)
	subj1, subj2 = tp1[0], tp2[0]
	order1, order2 = int(tp1[1]), int(tp2[1])
	if cohort == 'backhed':
		ilabels = ['I_Birth', 'I_4mon','I_12mon']
		mlabels = ['M_Birth']
	elif cohort == 'ferretti':
		ilabels = ['I_1day','I_3days', 'I_1wk', 'I_1mon', 'I_4mon']
		mlabels = ['M_Birth']
	elif cohort == 'yassour':
		ilabels = ['I_Birth', 'I_2wk','I_1mon', 'I_2mon', 'I_3mon']
		mlabels = ['M_Gest', 'M_Birth', 'M_3mon']
	else:
		return 'error: bad cohort'
	label1 = mlabels[order1-1] if subj1 == 'M' else ilabels[order1-1]
	label2 = mlabels[order2-1] if subj2 == 'M' else ilabels[order2-1]
	return (label1.replace('_','-') + " > " + label2.replace('_','-'))

# Survival curve computation
def calculate_unnormalized_survival_from_vector(counts):
	counts = sorted(counts)
	xs = []
	ns = []
	ns_cur = len(counts)
	min = -1
	for count in counts:
		if count > min:
			ns.append(ns_cur) # Number of elements greater or equal
			xs.append(count)
			min = count
		ns_cur -= 1
	xs.append(xs[len(xs)-1]+1)
	ns.append(0)
	return xs, numpy.array(ns)

# Not efficient but enough for this tiny data
def in_bin(number, range_tuple):
	return number >= range_tuple[0] and number <= range_tuple[1]

# Fractions instead of absolute counts (only use numeric arrays)
def normalize_array(arr):
	return numpy.array(arr) / float(numpy.sum(arr))

def order_tps(tp_list):
	ordered_tps = []
	order_map = {}
	for tp1, tp2 in tp_list:
		order1 = tp1[1]
		order2 = tp2[1]
		if tp1[0] == 'M':
			order1 = '0'
		if tp2[0] == 'M':
			order2 = '0'
		prev, next = (order1, order2) if order1 < order2 else (order2, order1)
		order_map[int(prev + next)] = frozenset((tp1, tp2))
	for key in sorted(order_map.keys()):
		ordered_tps.append(order_map[key])
	return ordered_tps

def tp_to_xrange(tp_pair):
	tp1, tp2 = tp_pair
	order1 = int(tp1[1])
	order2 = int(tp2[1])
	if tp1[0] == 'M':
		order1 = 0
	if tp2[0] == 'M':
		order2 = 0
	prev, next = (order1, order2) if order1 < order2 else (order2, order1)
	return (prev, next)

# Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--sweep_type", help="Full or partial sweep", default="full")
parser.add_argument("--force", help="Forces regeneration of pickled data", default=False)
args = parser.parse_args()

# Full (0.6) or partial (0.3) sweep: full, partial
force = args.force if args.force == False else True
sweep_type = args.sweep_type
if sweep_type not in ['full', 'partial']:
	sys.exit("Invalid sweep_type. Choose from full, partial")
if sweep_type == 'full':
	lower_threshold, upper_threshold = 0.2, 0.8
elif sweep_type == 'partial':
	lower_threshold, upper_threshold = 0.35, 0.65

# Parameters
mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']	= False
mpl.rcParams['legend.fontsize']	 = 'small'

min_sample_size = 3
variant_types = ['1D','4D']
within_host_type = 'consecutive' # consecutive timepoints
#within_host_type = 'longest' # longest interval available
min_snp_change_sample_size = 5

modification_difference_threshold = config.modification_difference_threshold
replacement_difference_threshold = config.replacement_difference_threshold

# Cohorts
cohorts = ['backhed', 'ferretti', 'yassour', 'hmp']
mi_cohorts = ['backhed', 'ferretti', 'yassour']

# Pickle
pooled_snp_pickle_fn = '%s/pickles/pooled_snp_change_%s.pkl' % (config.data_directory, sweep_type)
pooled_gene_pickle_fn = '%s/pickles/pooled_gene_change_%s.pkl' % (config.data_directory, sweep_type)

pooled_between_snp_pickle_fn = '%s/pickles/pooled_between_snp_change_%s.pkl' % (config.data_directory, sweep_type)
pooled_between_gene_pickle_fn = '%s/pickles/pooled_between_gene_change_%s.pkl' % (config.data_directory, sweep_type)


if force == False and path.exists(pooled_snp_pickle_fn) and path.exists(pooled_between_snp_pickle_fn):
	pooled_snp_change_distribution = pickle.load(open(pooled_snp_pickle_fn, 'rb'))
	pooled_gene_change_distribution = pickle.load(open(pooled_gene_pickle_fn, 'rb'))
	
	pooled_between_snp_change_distribution = pickle.load(open(pooled_between_snp_pickle_fn, 'rb'))
	pooled_between_gene_change_distribution = pickle.load(open(pooled_between_gene_pickle_fn, 'rb'))
else:
	# Load subject and sample metadata
	sys.stderr.write("Loading sample metadata...\n")
	subject_sample_map = su.parse_subject_sample_map()
	sample_order_map = su.parse_sample_order_map()
	sample_country_map = su.parse_sample_country_map()
	sample_subject_map = su.calculate_sample_subject_map(subject_sample_map)
	all_samples = sample_country_map.keys()
	sys.stderr.write("Done!\n")
	
	# HMP samples
	hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

	# All mother-infant samples
	all_mi_samples = su.get_sample_names('all','all')
	mi_samples = {cohort: su.get_sample_names(cohort, 'all') for cohort in mi_cohorts}
	mother_samples = su.get_sample_names('mother','all')
	infant_samples = su.get_sample_names('infant','all')

	# Species list
	good_species_list = parse_midas_data.parse_good_species_list()

	# Species SNP/gene change distributions
	# Nested dictionaries: timepoint pair followed by species
	species_snp_change_distribution = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}
	species_snp_nerrs = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}
	species_gene_change_distribution = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}
	species_gene_nerrs = {cohort: defaultdict(lambda: defaultdict(list)) for cohort in cohorts}

	# Pooled SNP/gene change distributions
	pooled_snp_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
	pooled_gene_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
	pooled_snp_length_distribution = {cohort: defaultdict(list) for cohort in cohorts}
	pooled_gene_length_distribution = {cohort: defaultdict(list) for cohort in cohorts}

	# typical value, median other sample
	pooled_between_snp_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
	pooled_between_gene_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}

	# closest other sample
	pooled_min_between_snp_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}
	pooled_min_between_gene_change_distribution = {cohort: defaultdict(list) for cohort in cohorts}

	for species_name in good_species_list:
		
		sys.stderr.write("\nProcessing %s...\n" % species_name)
		
		qp_sample_sets = {}
		for cohort in cohorts:
			if cohort == 'hmp':
				qp_sample_sets[cohort] = su.calculate_qp_samples(hmp_samples, species_name)['qp']
			else:
				qp_sample_sets[cohort] = su.calculate_qp_samples(su.get_sample_names(cohort,'all'), species_name)['qp']
		qp_sample_lists = {cohort: sorted(qp_sample_sets[cohort]) for cohort in cohorts}
		combined_qp_samples = sorted(su.flatten([qp_sample_lists[cohort] for cohort in cohorts]))
		combined_sample_idx_map = {combined_qp_samples[i] : i for i in xrange(0,len(combined_qp_samples))}
		
		# Using Backhed to threshold on sample size
		sample_size = len(qp_sample_sets['backhed'])
		if sample_size < min_sample_size:
			continue
		
		import calculate_private_snvs
		private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)
		
		import calculate_snp_prevalences
		snv_freq_map = calculate_snp_prevalences.parse_population_freqs("all", species_name,polarize_by_consensus=True)
		snv_freq_keys = snv_freq_map.keys()
		snv_freq_values = snv_freq_map.values()
		
		import core_gene_utils
		gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
		gene_freq_values = numpy.array(gene_freq_map.values())
		gene_freq_weights = gene_freq_values*1.0/gene_freq_values.sum()
		
		# Combined Backhed, Ferretti, Yassour calculations
		sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
		substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
		if substitution_rate_map == {}: # Not enough haploid samples
			sys.stderr.write("Not enough haploid samples!\n")
			continue
		sys.stderr.write("Calculating SNV matrix...\n")
		dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=combined_qp_samples)
		
		snp_difference_matrix = snp_mut_difference_matrix+snp_rev_difference_matrix
		snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
		snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading gene matrix...\n")
		gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=combined_qp_samples)
		gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
		gene_opportunity_matrix = gene_loss_opportunity_matrix
		gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading 1D & 4D opportunity matrices...\n")
		
		difference_matrices = {}
		opportunity_matrices = {}
		
		for var_type in variant_types:		
			dummy_samples, difference_matrix, opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, var_type, allowed_samples=combined_qp_samples)
			
			difference_matrices[var_type] = difference_matrix
			opportunity_matrices[var_type] = opportunity_matrix
		
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading pre-computed temporal changes...\n")
		temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
		sys.stderr.write("Done!\n")
		
		### Now loop over different cohorts
		for cohort in cohorts:		
			desired_samples = qp_sample_lists[cohort]
			
			same_subject_idxs = su.calculate_mi_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)
			
			'''
			same_sample_idxs, same_subject_idxs, diff_subject_idxs = su.calculate_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type)
			
			# Include same mother and infant comparisons
			more_same_subject_idxs_i = []
			more_same_subject_idxs_j = []
			
			for i in range(len(desired_samples)):
				for j in range(i+1, len(desired_samples)):
					sample_i = desired_samples[i]
					sample_j = desired_samples[j]
					if su.is_mi_pair(sample_i, sample_j, mother_samples, infant_samples) and su.is_same_mi_subject(sample_i, sample_j, sample_subject_map):
						more_same_subject_idxs_i.append(i)
						more_same_subject_idxs_j.append(j)
			
			idxs1 = numpy.array(numpy.append(same_subject_idxs[0], more_same_subject_idxs_i), dtype=numpy.int32)
			idxs2 = numpy.array(numpy.append(same_subject_idxs[1], more_same_subject_idxs_j), dtype=numpy.int32)
			same_subject_idxs = (idxs1, idxs2)
			
			#apply_sample_index_map_to_indices(sample_idx_map, idxs):
			#new_idxs = (numpy.array([sample_idx_map[i] for i in idxs[0]]), numpy.array([sample_idx_map[i] for i in idxs[1]]))
			'''
			
			# Loop over different pairs of within-host samples
			for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):			 
					
					sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
					sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
					
					if sample_i in hmp_samples:
						tp_i = 'A' + str(sample_order_map[sample_i][1])
					else:
						tp_i = ('M' if sample_i in mother_samples else 'I') + str(sample_order_map[sample_i][1])
					if sample_j in hmp_samples:
						tp_j = 'A' + str(sample_order_map[sample_j][1])
					else:
						tp_j = ('M' if sample_j in mother_samples else 'I') + str(sample_order_map[sample_j][1])
					tp_pair = frozenset((tp_i, tp_j))
					
					i = combined_sample_idx_map[sample_i]
					j = combined_sample_idx_map[sample_j]
					
					good_idxs = su.calculate_samples_in_different_subjects(subject_sample_map, combined_qp_samples, sample_i)
					good_idxs *= ( (snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5) )
					
					if good_idxs.sum() < 1:
							continue
					
					L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
					
					if L<config.min_opportunities:
							continue
					
					nerr = L*perr
					
					num_mutations = len(mutations)
					num_reversions = len(reversions)
					num_snp_changes = num_mutations+num_reversions
					
					gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
					
					gene_nerr = gene_L*gene_perr
					num_gains = len(gains)
					num_losses = len(losses)
					num_gene_changes = num_gains+num_losses
					
					if (perr<-0.5) or (gene_perr < -0.5):
							continue
					
					if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
							continue # Only take things with low-ish FPR
					
					# Species specific distributions
					species_snp_change_distribution[cohort][tp_pair][species_name].append( num_snp_changes)
					species_snp_nerrs[cohort][tp_pair][species_name].append(nerr)
					
					species_gene_change_distribution[cohort][tp_pair][species_name].append(num_gene_changes)
					species_gene_nerrs[cohort][tp_pair][species_name].append(gene_nerr)
					
					# Pooled distributions
					pooled_snp_change_distribution[cohort][tp_pair].append(num_snp_changes)
					pooled_gene_change_distribution[cohort][tp_pair].append(num_gene_changes)
					
					pooled_snp_length_distribution[cohort][tp_pair].append(L)
					pooled_gene_length_distribution[cohort][tp_pair].append(gene_L)
					
					# Matched between-host samples
					# typical
					pooled_between_snp_change_distribution[cohort][tp_pair].append( choice(snp_difference_matrix[i, good_idxs]) ) 
					pooled_between_gene_change_distribution[cohort][tp_pair].append( choice(gene_difference_matrix[i, good_idxs]) )
					
					# minimum
					pooled_min_between_snp_change_distribution[cohort][tp_pair].append( snp_difference_matrix[i, good_idxs].min() )
					pooled_min_between_gene_change_distribution[cohort][tp_pair].append( gene_difference_matrix[i, good_idxs].min() )
	
	pickle.dump(pooled_snp_change_distribution, open(pooled_snp_pickle_fn, 'w'))
	pickle.dump(pooled_gene_change_distribution, open(pooled_gene_pickle_fn, 'w'))
	pickle.dump(pooled_between_snp_change_distribution, open(pooled_between_snp_pickle_fn, 'w'))
	pickle.dump(pooled_between_gene_change_distribution, open(pooled_between_gene_pickle_fn, 'w'))

###################################
#
# Set up figure 1: Distribution of SNV changes
#
###################################

fig_snp, ax_snp = plt.subplots(1, 3, sharey='row', figsize=(30,6))

colormap = cmx.get_cmap('jet', 10)
colors = [colormap(x) for x in numpy.array([x for x in range(0,10)])/10.0]

custom_backhed_tp_pairs = [frozenset(['I1', 'M1']),  frozenset(['I2', 'M1']), frozenset(['I3', 'M1']), frozenset(['I1', 'I2']), frozenset(['I3', 'I2']), frozenset(['I1', 'I3'])]

# Plot SNP change distribution

for i in range(len(mi_cohorts)):
	cohort = mi_cohorts[i]
	ax_snp[i].set_xscale('log')
	ax_snp[i].set_yscale('log')
	ax_snp[i].set_ylabel('Fraction comparisons $\geq n$', fontsize=11)
	ax_snp[i].set_xlabel('# SNP changes', fontsize=11)
	
	ax_snp[i].spines['top'].set_visible(False)
	ax_snp[i].spines['right'].set_visible(False)
	ax_snp[i].get_xaxis().tick_bottom()
	ax_snp[i].get_yaxis().tick_left()
	
	color_i = 0
	if cohort == 'backhed':
		# Within-host, adult
		counts = []
		for tp_pair in pooled_snp_change_distribution['hmp'].keys():
			counts += pooled_snp_change_distribution['hmp'][tp_pair]
		xs, ns = calculate_unnormalized_survival_from_vector(counts)
		mlabel = 'HMP: Adult 6 months' + (' (n=%d)' % ns[0])
		ax_snp[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
		ymin = 1.0/ns[0]
		ymax = 1.3
		ax_snp[i].set_ylim([ymin,ymax])
		color_i += 1
		# Unrelated adults
		counts = []
		for tp_pair in pooled_between_snp_change_distribution['hmp'].keys():
			counts += pooled_between_snp_change_distribution['hmp'][tp_pair]
		xs, ns = calculate_unnormalized_survival_from_vector(counts)
		ax_snp[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label="Unrelated adults", where='pre',zorder=4)
		color_i += 1
		# Within-host, mi
		for tp_pair in custom_backhed_tp_pairs:
			counts = pooled_snp_change_distribution[cohort][tp_pair]
			if len(counts) < min_snp_change_sample_size:
				continue		
			xs, ns = calculate_unnormalized_survival_from_vector(counts)
			mlabel = format_tp_label(cohort, tp_pair) + (' (n=%d)' % ns[0])
			ax_snp[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
			color_i += 1
	else:
		# Within-host
		for tp_pair in pooled_snp_change_distribution[cohort].keys():
			counts = pooled_snp_change_distribution[cohort][tp_pair]		
			if len(counts) < min_snp_change_sample_size:
				continue		
			xs, ns = calculate_unnormalized_survival_from_vector(counts)
			mlabel = format_tp_label(cohort, tp_pair) + (' (n=%d)' % ns[0])
			ax_snp[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
			color_i += 1	
	# Unrelated
	counts = []
	for tp_pair in pooled_between_snp_change_distribution[cohort].keys():
		counts += pooled_between_snp_change_distribution[cohort][tp_pair]
	xs, ns = calculate_unnormalized_survival_from_vector(counts)
	ax_snp[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label="Unrelated mother/infants", where='pre',zorder=4)
	
	ax_snp[i].legend(loc='best', frameon=True, fontsize=10, numpoints=1, ncol=1, handlelength=1)
	
	if sweep_type == "partial":
		modification_difference_threshold = 50
	
	# Now fill in the graphics
	ax_snp[i].fill_between([1e-01,1], [ymin,ymin],[ymax,ymax],color='0.8',zorder=1)
	ax_snp[i].fill_between([1e0,modification_difference_threshold],[ymin,ymin],[ymax,ymax],color='#deebf7',zorder=1)
	ax_snp[i].fill_between([replacement_difference_threshold,1e05],[ymin,ymin],[ymax,ymax],color='#fee0d2',zorder=1)
	
	ax_snp[i].text( exp((log(1e05)+log(replacement_difference_threshold))/2), ymax*1.2, 'putative\nreplacement',fontsize=12,fontstyle='italic',ha='center',color='#fc9272',zorder=1)
	ax_snp[i].text( exp((log(1)+log(modification_difference_threshold))/2), ymax*1.2, 'putative\nmodification',fontsize=12,fontstyle='italic',ha='center',color='#9ecae1',zorder=1)

###################################
#
# Only plot Backhed
#
###################################

fig_b_snp, ax_b_snp = plt.subplots(figsize=(12,6))
colors = ['#9c1642','#2a6626','#74ad3b','#a9db1f','#0300ad','#145dc9','#0cdff2']

ax_b_snp.set_xscale('log')
ax_b_snp.set_yscale('log')
ax_b_snp.set_ylabel('Fraction comparisons $\geq n$', fontsize=11)
ax_b_snp.set_xlabel('# SNP changes', fontsize=11)
ax_b_snp.spines['top'].set_visible(False)
ax_b_snp.spines['right'].set_visible(False)
ax_b_snp.get_xaxis().tick_bottom()
ax_b_snp.get_yaxis().tick_left()

color_i = 0
counts = []
for tp_pair in pooled_snp_change_distribution['hmp'].keys():
	counts += pooled_snp_change_distribution['hmp'][tp_pair]
xs, ns = calculate_unnormalized_survival_from_vector(counts)
mlabel = 'HMP: Adult 6mon' + ('\n(n=%d)' % ns[0])
ax_b_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
ymin = 1.0/ns[0]
ymax = 1.3
ax_b_snp.set_ylim([ymin,ymax])
color_i += 1
for tp_pair in custom_backhed_tp_pairs:
	counts = pooled_snp_change_distribution['backhed'][tp_pair]		
	if len(counts) < min_snp_change_sample_size:
		continue		
	xs, ns = calculate_unnormalized_survival_from_vector(counts)
	mlabel = format_tp_label('backhed', tp_pair) + ('\n(n=%d)' % ns[0])
	ax_b_snp.step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
	color_i += 1

chartBox = ax_b_snp.get_position()
ax_b_snp.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.8, chartBox.height])
leg = ax_b_snp.legend(loc='upper center', bbox_to_anchor=(1.2, 0.8), frameon=False, fontsize=12, numpoints=1, ncol=1, handlelength=1)
for legobj in leg.legendHandles:
	legobj.set_linewidth(2.2)

# Now fill in the graphics
ax_b_snp.fill_between([1e-01,1], [ymin,ymin],[ymax,ymax],color='0.8',zorder=1)
ax_b_snp.fill_between([1e0,modification_difference_threshold],[ymin,ymin],[ymax,ymax],color='#deebf7',zorder=1)
ax_b_snp.fill_between([replacement_difference_threshold,1e05],[ymin,ymin],[ymax,ymax],color='#fee0d2',zorder=1)
ax_b_snp.text( exp((log(1e05)+log(replacement_difference_threshold))/2), ymax*1.2, 'putative\nreplacement',fontsize=12,fontstyle='italic',ha='center',color='#fc9272',zorder=1)
ax_b_snp.text( exp((log(1)+log(modification_difference_threshold))/2), ymax*1.2, 'putative\nmodification',fontsize=12,fontstyle='italic',ha='center',color='#9ecae1',zorder=1)

###################################
#
# Set up figure 2: Distribution of gene changes
#
###################################

fig_gene, ax_gene = plt.subplots(1, 3, sharey='row', figsize=(24,6))

colormap = cmx.get_cmap('jet', 10)
colors = [colormap(x) for x in numpy.array([x for x in range(0,10)])/10.0]

# Plot gene change distribution

for i in range(len(mi_cohorts)):
	cohort = mi_cohorts[i]
	ax_gene[i].loglog([0.1],[1],'k.')
	ax_gene[i].set_ylabel('Fraction comparisons $\geq n$')
	ax_gene[i].set_xlabel('# gene changes')
	# ax_gene[i].set_xlim([0.6,1e05])
	
	ax_gene[i].spines['top'].set_visible(False)
	ax_gene[i].spines['right'].set_visible(False)
	ax_gene[i].get_xaxis().tick_bottom()
	ax_gene[i].get_yaxis().tick_left()
	
	color_i = 0
	for tp_pair in pooled_gene_change_distribution[cohort].keys():
		counts = pooled_gene_change_distribution[cohort][tp_pair]		
		if len(counts) < min_snp_change_sample_size:
			continue		
		xs, ns = calculate_unnormalized_survival_from_vector(counts)		
		mlabel = format_tp_label(cohort, tp_pair) + (' (n=%d)' % ns[0])
		ax_gene[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1, label=mlabel, where='pre',zorder=4)
		ymin = 1.0/ns[0]
		ymax = 1.3
		ax_gene[i].set_ylim([ymin,ymax])
		color_i += 1
	if cohort == 'backhed':
		counts = []
		for tp_pair in pooled_gene_change_distribution['hmp'].keys():
			counts += pooled_gene_change_distribution['hmp'][tp_pair]
		xs, ns = calculate_unnormalized_survival_from_vector(counts)
		mlabel = 'HMP: Adult 6 months' + (' (n=%d)' % ns[0])
		ax_gene[i].step(xs,ns/float(ns[0]),'-',color=colors[color_i],linewidth=1.4, label=mlabel, where='pre',zorder=4)
		ymin = 1.0/ns[0]
		ymax = 1.3
		ax_gene[i].set_ylim([ymin,ymax])
		color_i += 1
	
	ax_gene[i].legend(loc='upper right',frameon=False,fontsize=10,numpoints=1,ncol=1,handlelength=1)

###################################
#
# Set up figure 3: Distribution of SNV changes, alternate bar visualization
#
###################################
'''
fig_snp_bar = {} # Dictionary of each cohort's figure

colormap = cmx.get_cmap('spectral', 26)
colors = [colormap(x) for x in numpy.array([x for x in range(0,26)])/26.0]

bins = [(0,0),(1,3),(4,10),(11,99),(100,999),(1000,9999),(10000,9999999)]
binlabels = ['0', '1-3', '4-10', '11-99', '100-999', '1,000-9,999', '>10,000']

for cohort in mi_cohorts:
	snps_binned = {}
	x_pos = []
	tp_to_show = pooled_snp_change_distribution[cohort].keys()
	num_tp = len(tp_to_show)
	barWidth = 1.0/(num_tp + 2)
	
	fig, ax = plt.subplots(figsize=(num_tp*2,4))
	fig_snp_bar[cohort] = fig
	
	for i in range(num_tp):
		tp_pair = tp_to_show[i]
		snps_binned[tp_pair] = defaultdict(int)
		for snv_change in pooled_snp_change_distribution[cohort][tp_pair]:
			for bin in bins:
				if in_bin(snv_change, bin):
					snps_binned[tp_pair][bin] += 1
		if i == 0:
			x_pos.append(barWidth + numpy.arange(len(bins)))
		else:
			x_pos.append([x + barWidth for x in x_pos[i-1]])
	
	for i in range(len(bins)):
		plt.axvline(x=i, color='#9c9c9c')
	
	for i in range(num_tp):
		tp_pair = tp_to_show[i]
		list = []
		for bin in bins:
			list.append(snps_binned[tp_pair][bin])
		mlabel = format_tp_label(cohort, tp_pair) + (' (n=%d)' % len(pooled_snp_change_distribution[cohort][tp_pair]))
		ax.bar(x_pos[i], list, color=colors[i], width=barWidth, edgecolor='white', label=mlabel)
	
	ax.set_xlabel('Number of SNP changes', fontweight='bold')
	ax.set_ylabel('Number of sample pairs', fontweight='bold')
	ax.set_xticks(0.5 + numpy.arange(len(bins)))
	ax.set_xticklabels(binlabels)
	ax.legend(loc='best', title='Timepoint pair',fontsize='medium')
	ax.set_title(cohort[0].upper() + cohort[1:])
'''
###################################
#
# Set up figure 4: Distribution of SNV changes, alt-alt bar visualization
#
###################################
'''
fig_snp_bar_alt = {} # Dictionary of each cohort's figure, bar plot

colormap = cmx.get_cmap('coolwarm', 10)
colors = [colormap(x) for x in numpy.array([x for x in range(0,7)])/10.0]

bins = [(0,0),(1,3),(4,10),(11,99),(100,999),(1000,9999),(10000,9999999)]
binlabels = ['0', '1-3', '4-10', '11-99', '100-999', '1,000-9,999', '>10,000']

for cohort in mi_cohorts:
	snps_binned = {}
	x_pos = []
	if cohort == 'backhed':
		tp_to_show = [frozenset(['M1', 'I1']), frozenset(['I1', 'I2']), frozenset(['I2', 'I3'])]
	elif cohort == 'ferretti':
		tp_to_show = [frozenset(['M1', 'I1']), frozenset(['I1', 'I2']), frozenset(['I2', 'I3']), frozenset(['I3', 'I4']), frozenset(['I4', 'I5'])]
	elif cohort == 'yassour':
		tp_to_show = [frozenset(['M2', 'I1']), frozenset(['I1', 'I2']), frozenset(['I2', 'I3']), frozenset(['I3', 'I4']), frozenset(['I4', 'I5'])]
	num_tp = len(tp_to_show)
	barWidth = 1.0/(len(bins) + 2)
	
	tp_labels = [format_tp_label(cohort, tp_pair) + (' (n=%d)' % len(pooled_snp_change_distribution[cohort][tp_pair])) for tp_pair in tp_to_show]
	
	fig, ax = plt.subplots(figsize=(num_tp*2,4))
	fig_snp_bar_alt[cohort] = fig
	
	for i in range(num_tp):
		tp_pair = tp_to_show[i]
		snps_binned[tp_pair] = defaultdict(int)
		for snv_change in pooled_snp_change_distribution[cohort][tp_pair]:
			for bin in bins:
				if in_bin(snv_change, bin):
					snps_binned[tp_pair][bin] += 1
	
	for i in range(len(bins)):
		if i == 0:
			x_pos.append(barWidth + numpy.arange(num_tp))
		else:
			x_pos.append([x + barWidth for x in x_pos[i-1]])
	
	for i in range(len(bins)):
		bin = bins[i]
		list = []
		for tp_pair in tp_to_show:
			total_num = len(pooled_snp_change_distribution[cohort][tp_pair])
			if total_num != 0:
				list.append(snps_binned[tp_pair][bin]/float(total_num))
		ax.bar(x_pos[i], list, color=colors[i], width=barWidth, edgecolor='white', label=binlabels[i])
	
	for i in range(num_tp):
		plt.axvline(x=i, color='#9c9c9c')
	
	ax.set_xlabel('Timepoint pair', fontweight='bold')
	ax.set_ylabel('Proportion of sample pairs', fontweight='bold')
	ax.set_xticks(0.5 + numpy.arange(num_tp))
	ax.set_xticklabels(tp_labels)
	ax.legend(loc='best', title='Number of SNP changes',fontsize='medium')
	ax.set_title(cohort[0].upper() + cohort[1:])
'''
###################################
#
# Set up figure 5: Distribution of SNV changes, alt pie visualization
#
###################################
'''
fig_snp_pie = {}

for cohort in mi_cohorts:
	snps_binned = {}
	num_tp = len(pooled_snp_change_distribution[cohort].keys())

	tp_labels = [format_tp_label(cohort, tp_pair) + (' (n=%d)' % len(pooled_snp_change_distribution[cohort][tp_pair])) for tp_pair in pooled_snp_change_distribution[cohort].keys()]

	fig, ax = plt.subplots(1, num_tp, figsize=(num_tp*3,4))
	fig_snp_pie[cohort] = fig

	for i in range(num_tp):
		tp_pair = pooled_snp_change_distribution[cohort].keys()[i]
		snps_binned[tp_pair] = defaultdict(int)
		for snv_change in pooled_snp_change_distribution[cohort][tp_pair]:
			for bin in bins:
				if in_bin(snv_change, bin):
					snps_binned[tp_pair][bin] += 1

	for i in range(num_tp):
		tp_pair = pooled_snp_change_distribution[cohort].keys()[i]
		list = []
		for bin in bins:
			list.append(snps_binned[tp_pair][bin])
		ax[i].pie(list, colors=colors, labels=binlabels)
		ax[i].set_title(format_tp_label(cohort, tp_pair))
	
	fig.suptitle(cohort[0].upper() + cohort[1:])
'''
###################################
#
# Set up figure 6: Distribution of SNV changes, timeline visualization
#
###################################

fig_snp_tl = {}

colormap = cmx.get_cmap('coolwarm', 7)
colors = [colormap(x) for x in numpy.array([x for x in range(0,7)])/7.0]

bins = [(0,0),(1,3),(4,10),(11,99),(100,999),(1000,9999),(10000,9999999)]
bin_to_color = {bins[i]: colors[i] for i in range(len(bins))}
binlabels = ['0', '1-3', '4-10', '11-99', '100-999', '1,000-9,999', '>10,000']

timeline_labels = {'backhed': ['Mother: Delivery', '2-5 Days', '4 Months','12 Months'], 'ferretti': ['Mother: Delivery', '1 Day','3 Days', '1 Week', '1 Month', '4 Months'], 'yassour': ['Mother: Delivery', 'Birth', '2 Weeks','1 Month', '2 Months', '3 Months']}

cohort_spacing = {'backhed': 0.03, 'ferretti': 0.1, 'yassour': 0.03}

for cohort in mi_cohorts:
	snps_binned = {}
	tp_to_show = pooled_snp_change_distribution[cohort].keys()
	tp_to_show = [tp for tp in tp_to_show if ('M2' not in tp and 'M3' not in tp)]
	tp_to_show = order_tps(tp_to_show)
	tp_labels = [format_tp_label(cohort, tp_pair) + (' (n=%d)' % len(pooled_snp_change_distribution[cohort][tp_pair])) for tp_pair in tp_to_show]
	num_tp = len(tp_to_show)
	
	for tp_pair in tp_to_show:
		snps_binned[tp_pair] = defaultdict(int)
		for snv_change in pooled_snp_change_distribution[cohort][tp_pair]:
			for bin in bins:
				if in_bin(snv_change, bin):
					snps_binned[tp_pair][bin] += 1
	
	ns = [len(pooled_snp_change_distribution[cohort][tp_pair]) for tp_pair in tp_to_show]
	
	# fig, ax = plt.subplots(figsize=(7, sum(ns)*0.03))	
	spacing = cohort_spacing[cohort]
	fig, ax = plt.subplots(num_tp, 1, sharex = True, figsize=(7, sum(ns)*spacing), gridspec_kw={'height_ratios': ns, 'hspace': 0})
	fig_snp_tl[cohort] = fig
	
	ax[0].set_xlim(0, len(timeline_labels[cohort])-1)
	ax[0].set_xticks(range(len(timeline_labels[cohort])))
	ax[0].set_xticklabels(timeline_labels[cohort])
	ax[0].xaxis.set_label_position('top') 
	ax[0].xaxis.tick_top()
	
	for i in range(num_tp):
		tp_pair = tp_to_show[i]
		yrange = range(ns[i]-1, -1, -1)
		xmin, xmax = tp_to_xrange(tp_pair)
		color_list = []
		for bin in bins:
			for _ in range(snps_binned[tp_pair][bin]):
				color_list.append(bin_to_color[bin])
		lines = ax[i].hlines(yrange, xmin, xmax, colors=color_list)
		lines.set_linewidth(spacing*64)
		ax[i].set_ylim(0, ns[i])
		if ns[i] > 2:
			ax[i].set_ylabel('n = %s' % ns[i], rotation=0)
		ax[i].yaxis.set_label_coords(-0.05, 0.5)
		ax[i].yaxis.set_ticks_position('none')
		ax[i].set_yticklabels([])
		ax[i].spines['right'].set_visible(False)
		ax[i].spines['bottom'].set_visible(False)
		if i != 0:
			ax[i].spines['top'].set_visible(False)
			ax[i].xaxis.set_ticks_position('none')
	
	# Last axis
	legend_elements = [Patch(facecolor=colors[k], label=binlabels[k]) for k in range(len(bins))]
	ax[i].legend(handles=legend_elements, loc='lower left', title='Number of SNP changes', frameon=True)
	ax[0].set_title(cohort[0].upper() + cohort[1:], y = 1.2)

# Save figures

sys.stderr.write("Saving figures...\t")
fig_snp.savefig('%s/temporal_snp_changes_%s_pooled.pdf' % (config.analysis_directory, sweep_type),bbox_inches='tight')
fig_b_snp.savefig('%s/temporal_snp_changes_%s_pooled_backhed.pdf' % (config.analysis_directory, sweep_type),bbox_inches='tight')
fig_gene.savefig('%s/temporal_gene_changes_%s_pooled.pdf' % (config.analysis_directory, sweep_type),bbox_inches='tight')
for cohort in mi_cohorts:
	# fig_snp_bar[cohort].savefig('%s/temporal_snp_changes_%s_pooled_bar_%s.pdf' % (config.analysis_directory, sweep_type, cohort),bbox_inches='tight')
	# fig_snp_bar_alt[cohort].savefig('%s/temporal_snp_changes_%s_pooled_bar_alt_clean_%s.pdf' % (config.analysis_directory, sweep_type, cohort),bbox_inches='tight')
	# fig_snp_pie[cohort].savefig('%s/temporal_snp_changes_%s_pooled_pie_%s.pdf' % (config.analysis_directory, sweep_type, cohort),bbox_inches='tight')
	fig_snp_tl[cohort].savefig('%s/temporal_snp_changes_%s_pooled_timeline_%s.pdf' % (config.analysis_directory, sweep_type, cohort),bbox_inches='tight')
sys.stderr.write("Done!\n")
