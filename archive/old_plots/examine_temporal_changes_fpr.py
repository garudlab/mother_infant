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
import calculate_singletons
import sfs_utils
from collections import defaultdict
from numpy.random import choice, binomial

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
	elif cohort == 'hmp':
		return (str(order1) + " > " + str(order2))
	else:
		return 'error: bad cohort'
	label1 = mlabels[order1-1] if subj1 == 'M' else ilabels[order1-1]
	label2 = mlabels[order2-1] if subj2 == 'M' else ilabels[order2-1]
	return (label1.replace('_','-') + " > " + label2.replace('_','-'))

# Convert timepoint pair frozenset to pair type name
def timepoint_pair_type(tp_pair):
	tp_pair_name_map = {'mother-mother': 'MM', 'infant-infant': 'II', 'mother-infant': 'MI', 'adult-adult': 'AA'}
	
	tp1, tp2 = tp_pair
	s1, s2 = tp1[0], tp2[0]
	
	if s1 == s2:
		return (s1 + s2) # MM, II or AA
	elif (s1 + s2) == 'MI' or (s1 + s2) == 'IM':
		return 'MI'
	else: return 'Invalid timepoint pair'

# Cohorts
cohorts = ['backhed', 'ferretti', 'yassour', 'hmp']
mi_cohorts = ['backhed', 'ferretti', 'yassour']

# Timepoint pair types
tp_pair_names = ['MM', 'MI', 'II', 'AA']

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.calculate_sample_subject_map(subject_sample_map)
sys.stderr.write("Done!\n")

# Pickled QP data
all_samples_pickle_fn = '%s/pickles/all_samples.pkl' % config.data_directory
all_samples = pickle.load(open(all_samples_pickle_fn, 'rb'))

qp_samples_pickle_fn = '%s/pickles/qp_samples.pkl' % config.data_directory
qp_samples = pickle.load(open(qp_samples_pickle_fn, 'rb'))

# HMP samples
hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

# All mother-infant samples
all_mi_samples = su.get_sample_names('all','all')
mi_samples = {cohort: su.get_sample_names(cohort, 'all') for cohort in mi_cohorts}
mother_samples = su.get_sample_names('mother','all')
infant_samples = su.get_sample_names('infant','all')

# Species list
good_species_list = parse_midas_data.parse_good_species_list()

# Parameters
sweep_types = {'full': (0.2, 0.8), 'partial': (0.35, 0.65)}
min_sample_size = 3

# Dictionaries containing cohort -> sweep type -> species -> tp_pair -> perr, nerr
# Include pooled in addition to species ('all_species') and in addition to tp_pair ('all_tp')
snps_fpr = {cohort: {stype: {} for stype in sweep_types.keys()} for cohort in cohorts}
genes_fpr = {cohort: {stype: {} for stype in sweep_types.keys()} for cohort in cohorts}

# Dictionaries containing cohort -> gene name -> list of perr/nerr
perr_all_snp_changes = {cohort: defaultdict(list) for cohort in cohorts}
nerr_all_snp_changes = {cohort: defaultdict(list) for cohort in cohorts}

perr_all_gene_gainloss = {cohort: defaultdict(list) for cohort in cohorts}
nerr_all_gene_gainloss = {cohort: defaultdict(list) for cohort in cohorts}

for species_name in good_species_list:
		sys.stderr.write("\nProcessing %s...\n" % species_name)	
		qp_sample_sets = {cohort: [all_samples[index] for index in qp_samples[cohort][species_name]] for cohort in cohorts}
		qp_sample_lists = {cohort: sorted(qp_sample_sets[cohort]) for cohort in cohorts}
		combined_qp_samples = sorted(su.flatten([qp_sample_lists[cohort] for cohort in cohorts]))
		combined_sample_idx_map = {combined_qp_samples[i] : i for i in xrange(0,len(combined_qp_samples))}
		
		# Using Backhed to threshold on sample size
		sample_size = len(qp_sample_sets['backhed'])
		if sample_size < min_sample_size:
			continue
		
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
		
		sys.stderr.write("Loading pre-computed temporal changes...\n")
		temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
		sys.stderr.write("Done!\n")
		
		### Now loop over different cohorts
		for cohort in cohorts:		
			desired_samples = qp_sample_lists[cohort]
			
			within_host_type_param = 'nonconsecutive'
			
			same_subject_idxs = su.calculate_mi_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type_param, one_per_mi_pair=False)
			
			for stype in sweep_types.keys():
				snps_fpr[cohort][stype][species_name] = {}
				genes_fpr[cohort][stype][species_name] = {}
			
			# Loop over different pairs of within-host samples
			for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):			 
				
				sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]]
				sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
				
				subject = sample_order_map[sample_i][0][:-2]
				
				if sample_i and sample_j in hmp_samples:
					tp_i, tp_j = ('A' + str(sample_order_map[sample_i][1]), 'A' + str(sample_order_map[sample_j][1]))
				else:
					tp_i = ('M' if sample_i in mother_samples else 'I') + str(sample_order_map[sample_i][1])
					tp_j = ('M' if sample_j in mother_samples else 'I') + str(sample_order_map[sample_j][1])
				
				tp_pair = frozenset((tp_i, tp_j))
				tp_pair_type = timepoint_pair_type(tp_pair)
				
				i = combined_sample_idx_map[sample_i]
				j = combined_sample_idx_map[sample_j]
				
				# print(sample_order_map[sample_i][0] + " " + str(sample_order_map[sample_i][1]) + " > " + sample_order_map[sample_j][0] + " " + str(sample_order_map[sample_j][1]))
				
				good_idxs = su.calculate_samples_in_different_subjects(subject_sample_map, combined_qp_samples, sample_i)
				good_idxs *= ((snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5))
				
				if good_idxs.sum() < 1:
						continue
				
				err_measures = ['perr', 'nerr']
				
				# Record perr and nerr for snp changes
				for stype in sweep_types.keys():
					
					lower_threshold, upper_threshold = sweep_types[stype]
					
					L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
					
					# Only ignore if not enough opportunities
					if L<config.min_opportunities:
						continue
					
					nerr = L*perr
					
					# Save gene name - err information
					all_snp_changes = mutations+reversions					
					for snp_change in all_snp_changes:
							gene_name = snp_change[0]
							perr_all_snp_changes[cohort][gene_name].append(perr)
							nerr_all_snp_changes[cohort][gene_name].append(nerr)
					
					# Individual
					if tp_pair not in snps_fpr[cohort][stype][species_name]:
						snps_fpr[cohort][stype][species_name][tp_pair] = {meas: [] for meas in err_measures}
				
					snps_fpr[cohort][stype][species_name][tp_pair]['perr'].append(perr)
					snps_fpr[cohort][stype][species_name][tp_pair]['nerr'].append(nerr)
					
					# Pooled across species
					if 'all_species' not in snps_fpr[cohort][stype]:
						snps_fpr[cohort][stype]['all_species'] = {}
					
					if tp_pair not in snps_fpr[cohort][stype]['all_species']:
						snps_fpr[cohort][stype]['all_species'][tp_pair] = {meas: [] for meas in err_measures}
					
					snps_fpr[cohort][stype]['all_species'][tp_pair]['perr'].append(perr)
					snps_fpr[cohort][stype]['all_species'][tp_pair]['nerr'].append(nerr)
					
					# Pooled across timepoint pair types
					if 'all_tp' not in snps_fpr[cohort][stype][species_name]:
						snps_fpr[cohort][stype][species_name]['all_tp'] = {meas: [] for meas in err_measures}
					
					snps_fpr[cohort][stype][species_name]['all_tp']['perr'].append(perr)
					snps_fpr[cohort][stype][species_name]['all_tp']['nerr'].append(nerr)
				
				# Record perr and nerr for gene changes
				gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
				
				gene_nerr = gene_L*gene_perr
				
				# Save gene name - err information
				all_gene_gainloss = gains+losses					
				for gene_gl in all_gene_gainloss:
						perr_all_gene_gainloss[cohort][gene_gl].append(gene_perr)
						nerr_all_gene_gainloss[cohort][gene_gl].append(gene_nerr)
				
				if tp_pair not in genes_fpr[cohort][stype][species_name]:
					genes_fpr[cohort][stype][species_name][tp_pair] = {meas: [] for meas in err_measures}
				
				genes_fpr[cohort][stype][species_name][tp_pair]['perr'].append(gene_perr)
				genes_fpr[cohort][stype][species_name][tp_pair]['nerr'].append(gene_nerr)
				
				# These thresholds may have to be adjusted, so I don't filter anything above
				'''
				if (perr<-0.5) or (gene_perr < -0.5):
						continue
				
				if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
						continue # Only take things with low-ish FPR
				'''

# Dump snps_fpr and genes_fpr
snps_fpr_pkl_fn = '%s/pickles/snps_fpr.pkl' % config.data_directory
genes_fpr_pkl_fn = '%s/pickles/genes_fpr.pkl' % config.data_directory
perr_by_gene_snp_pkl_fn = '%s/pickles/perr_by_gene_snp_change.pkl' % config.data_directory
nerr_by_gene_snp_pkl_fn = '%s/pickles/nerr_by_gene_snp_change.pkl' % config.data_directory
perr_by_gene_gainloss_pkl_fn = '%s/pickles/perr_by_gene_gainloss.pkl' % config.data_directory
nerr_by_gene_gainloss_pkl_fn = '%s/pickles/nerr_by_gene_gainloss.pkl' % config.data_directory

pickle.dump(perr_all_snp_changes, open(perr_by_gene_snp_pkl_fn, 'w'))
pickle.dump(nerr_all_snp_changes, open(nerr_by_gene_snp_pkl_fn, 'w'))
pickle.dump(perr_all_gene_gainloss, open(perr_by_gene_gainloss_pkl_fn, 'w'))
pickle.dump(nerr_all_gene_gainloss, open(nerr_by_gene_gainloss_pkl_fn, 'w'))
# pickle.dump(snps_fpr, open(snps_fpr_pkl_fn, 'w'))
# pickle.dump(genes_fpr, open(genes_fpr_pkl_fn, 'w'))

snps_fpr = pickle.load(open(snps_fpr_pkl_fn, 'r'))
genes_fpr = pickle.load(open(genes_fpr_pkl_fn, 'r'))
perr_all_snp_changes = pickle.load(open(perr_by_gene_snp_pkl_fn, 'r'))
nerr_all_snp_changes = pickle.load(open(nerr_by_gene_snp_pkl_fn, 'r'))
perr_all_gene_gainloss = pickle.load(open(perr_by_gene_gainloss_pkl_fn, 'r'))
nerr_all_gene_gainloss = pickle.load(open(nerr_by_gene_gainloss_pkl_fn, 'r'))

# Plot histograms of perrs and nerrs for each tuple (sweep_type, cohort)
fig_nerr, ax_nerr = plt.subplots(2, 3, figsize=(18,8)); i = 0; j = 0

for sweep_type in sweep_types.keys():
	j = 0
	for cohort in cohorts:
		snps_fpr_lists = snps_fpr[cohort][sweep_type]['all_species']
		snps_fpr_flat = [fpr for tp in snps_fpr_lists for fpr in snps_fpr_lists[tp]['nerr']]
		ax_nerr[i][j].hist(snps_fpr_flat, bins=20)
		ax_nerr[i][j].set_xlabel("Nerr values")
		ax_nerr[i][j].set_title("Nerr for %s, %s sweeps" % (cohort, sweep_type))
		ax_nerr[i][j].set_xscale('log')
		j += 1
	i += 1

fig_nerr.savefig('%s/nerr_hist.pdf' % config.analysis_directory, bbox_inches='tight')

# perr
fig_perr, ax_perr = plt.subplots(2, 3, figsize=(18,8)); i = 0; j = 0

for sweep_type in sweep_types.keys():
	j = 0
	for cohort in cohorts:
		snps_fpr_lists = snps_fpr[cohort][sweep_type]['all_species']
		snps_fpr_flat = [fpr for tp in snps_fpr_lists for fpr in snps_fpr_lists[tp]['perr']]		
		ax_perr[i][j].hist(snps_fpr_flat, bins=20)
		ax_perr[i][j].set_xlabel("Perr values")
		ax_perr[i][j].set_title("Perr for %s, %s sweeps" % (cohort, sweep_type))
		ax_perr[i][j].set_xscale('log')
		j += 1
	i += 1

fig_perr.savefig('%s/perr_hist.pdf' % config.analysis_directory, bbox_inches='tight')
