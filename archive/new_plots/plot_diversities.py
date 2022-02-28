from utils import sample_utils as su, config, parse_midas_data, stats_utils, sfs_utils
import pylab, sys, numpy as np, random, math
from utils import temporal_changes_utils
from collections import defaultdict
import bz2

import pickle
adir = config.analysis_directory
ddir = config.data_directory
pdir = "%s/pickles" % config.data_directory

# Load species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# =================================================================
# Load QP, polymorphism, alpha diversity info
sample_species_qp_dict = pickle.load(open("%s/sample_species_qp_dict.pkl" % pdir, 'rb'))
num_qp_dict = pickle.load(open("%s/plot_qp/num_qp_dict_binned.pkl" % pdir, 'rb'))
num_non_dict = pickle.load(open("%s/plot_qp/num_non_dict_binned.pkl" % pdir, 'rb'))
num_lowcov_dict = pickle.load(open("%s/plot_qp/num_lowcov_dict_binned.pkl" % pdir, 'rb'))
sample_species_polymorphism_dict = pickle.load(open("%s/sample_species_polymorphism_dict.pkl" % (pdir), 'rb'))
alpha_div_dict = pickle.load(open("%s/alpha_div_dict.pkl" % pdir, 'rb'))
# =================================================================

# =================================================================
# Load polymorphism (common and rare pS) info

sample_species_common_pS_dict = defaultdict(dict)
sample_species_rare_pS_dict = defaultdict(dict)
sample_species_seg_pS_dict = defaultdict(dict)

type, type_desc = 'common',  '0.2-0.5'
within_total_sites_nonQP = pickle.load(open("%s/within_total_sites_%s_nonQP.pkl" % (pdir, type), 'rb'))
within_total_sites_QP = pickle.load(open("%s/within_total_sites_%s_QP.pkl" % (pdir, type), 'rb'))

for within_total_sites_dict in [within_total_sites_QP, within_total_sites_nonQP]:
	for sample in within_total_sites_dict:
		for species in within_total_sites_dict[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_dict[sample][species]
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			sample_species_common_pS_dict[sample][species] = (w4D*1.0 + 1.0)/(t4D + 1.0)

type, type_desc = 'rare', '0-0.05 (excl. 0)'
within_total_sites_nonQP = pickle.load(open("%s/within_total_sites_%s_nonQP.pkl" % (pdir, type), 'rb'))
within_total_sites_QP = pickle.load(open("%s/within_total_sites_%s_QP.pkl" % (pdir, type), 'rb'))

for within_total_sites_dict in [within_total_sites_QP, within_total_sites_nonQP]:
	for sample in within_total_sites_dict:
		for species in within_total_sites_dict[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_dict[sample][species]
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			sample_species_rare_pS_dict[sample][species] = (w4D*1.0 + 1.0)/(t4D + 1.0)

type, type_desc = 'seg', '>0'
within_total_sites_nonQP = pickle.load(open("%s/within_total_sites_%s_nonQP.pkl" % (pdir, type), 'rb'))
within_total_sites_QP = pickle.load(open("%s/within_total_sites_%s_QP.pkl" % (pdir, type), 'rb'))

for within_total_sites_dict in [within_total_sites_QP, within_total_sites_nonQP]:
	for sample in within_total_sites_dict:
		for species in within_total_sites_dict[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_dict[sample][species]
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			sample_species_seg_pS_dict[sample][species] = (w4D*1.0 + 1.0)/(t4D + 1.0)
			if sample_species_seg_pS_dict[sample][species] > 0.7:
				print(sample + ", " + species)
				print(within_total_sites_dict[sample][species])
# =================================================================

# QP: aggregate over species for each timepoint

mi_tp_sample_dict, infant_tps_ordered = su.get_mi_tp_sample_dict(exclude_cohorts = ['olm'], binned = True)
mother_tps_ordered = sorted(mi_tp_sample_dict['mother'].keys())
tps_ordered_dict = {'mother': mother_tps_ordered, 'infant': infant_tps_ordered}

num_qp_agg_species = {'infant': [], 'mother': []}
num_non_agg_species = {'infant': [], 'mother': []}
num_lowcov_agg_species = {'infant': [], 'mother': []}

for cat in ['mother', 'infant']:
	for tp in tps_ordered_dict[cat]:
		
		if len(mi_tp_sample_dict[cat][tp]) < 10:
			continue # Skip timepoints with not enough data
		
		total_num_qp, total_num_non, total_num_lowcov = 0, 0, 0
	
		for species in good_species_list:
			total_num_qp += num_qp_dict[cat][species][tp]
			total_num_non += num_non_dict[cat][species][tp]
			total_num_lowcov += num_lowcov_dict[cat][species][tp]
		
		num_qp_agg_species[cat].append(total_num_qp)
		num_non_agg_species[cat].append(total_num_non)
		num_lowcov_agg_species[cat].append(total_num_lowcov)

# QP: Get proportion QP for each sample

num_qp_sample = {'infant': defaultdict(int), 'mother': defaultdict(int)}
num_nonqp_sample = {'infant': defaultdict(int), 'mother': defaultdict(int)}

for cat in ['mother', 'infant']:
	for sample in sample_species_qp_dict[cat]:
		for species in sample_species_qp_dict[cat][sample]:
			qp_status = sample_species_qp_dict[cat][sample][species]
			if qp_status == 'qp':
				num_qp_sample[cat][sample] += 1
			if qp_status == 'non-qp':
				num_nonqp_sample[cat][sample] += 1

prop_nonqp_sample = {'infant': [], 'mother': []}

for cat in ['mother', 'infant']:	
	for tp in tps_ordered_dict[cat]:
		
		num_samples = len(mi_tp_sample_dict[cat][tp])
		if num_samples < 10:
			continue # Skip timepoints with not enough data
		
		prop_nonqp_sample_tp = []
		
		for sample in mi_tp_sample_dict[cat][tp]:
			if sample in num_qp_sample[cat]:
				num_qp = num_qp_sample[cat][sample]
				num_nonqp = num_nonqp_sample[cat][sample]
				num_highcov = num_qp + num_nonqp
				if num_highcov >= 3:
					prop_nonqp = float(num_nonqp) / num_highcov
					prop_nonqp_sample_tp.append(prop_nonqp)		
				else:
					print("%s: Not enough high coverage species!" % sample)
		
		prop_nonqp_sample[cat].append(prop_nonqp_sample_tp)

# Plot time!

import numpy as np
from matplotlib import pyplot as plt

# =======================================================================
# Alpha diversity boxplot (all non-Olm infants) on top
# QP proportion on bottom
# =======================================================================

alpha_divs = [] # list of sample values for each tp
common_pSs = [] # list of sample values for each tp
rare_pSs = [] # list of sample values for each tp
seg_pSs = [] # list of sample values for each tp

labels = ['M:-3m', 'M:dlv', 'M:1d', 'M:2d', 'M:3m']

for i in range(len(mother_tps_ordered)):
	
	tp = mother_tps_ordered[i]
	num_samples = len(mi_tp_sample_dict['mother'][tp])
	if num_samples < 10:
		continue # Skip timepoints with not enough data
	
	labels[i] += ("\nn=%i" % num_samples)
	
	alpha_divs_tp = []
	common_pSs_tp	= []
	rare_pSs_tp = []
	seg_pSs_tp = []
	
	for sample in mi_tp_sample_dict['mother'][tp]:
		alpha_divs_tp.append(alpha_div_dict[sample])
		common_pSs_tp += [sample_species_common_pS_dict[sample][species] for species in sample_species_common_pS_dict[sample]]
		rare_pSs_tp += [sample_species_rare_pS_dict[sample][species] for species in sample_species_rare_pS_dict[sample]]
		seg_pSs_tp += [sample_species_seg_pS_dict[sample][species] for species in sample_species_seg_pS_dict[sample]]
	
	alpha_divs.append(alpha_divs_tp)
	common_pSs.append(common_pSs_tp)
	rare_pSs.append(rare_pSs_tp)
	seg_pSs.append(seg_pSs_tp)

for tp in infant_tps_ordered:
	
	num_samples = len(mi_tp_sample_dict['infant'][tp])
	if num_samples < 10:
		continue # Skip timepoints with not enough data
	
	labels.append(tp + "\n" + ("n=%i" % num_samples))
	
	alpha_divs_tp = []
	common_pSs_tp	= []
	rare_pSs_tp = []
	seg_pSs_tp = []
	
	for sample in mi_tp_sample_dict['infant'][tp]:
		alpha_divs_tp.append(alpha_div_dict[sample])
		common_pSs_tp += [sample_species_common_pS_dict[sample][species] for species in sample_species_common_pS_dict[sample]]
		rare_pSs_tp += [sample_species_rare_pS_dict[sample][species] for species in sample_species_rare_pS_dict[sample]]
		seg_pSs_tp += [sample_species_seg_pS_dict[sample][species] for species in sample_species_seg_pS_dict[sample]]
	
	alpha_divs.append(alpha_divs_tp)
	common_pSs.append(common_pSs_tp)
	rare_pSs.append(rare_pSs_tp)
	seg_pSs.append(seg_pSs_tp)

num_qp_infant = np.array(num_qp_agg_species['infant'])
num_non_infant = np.array(num_non_agg_species['infant'])
qp_props_infant = num_qp_infant/(num_non_infant + num_qp_infant).astype('float') # one value for each tp

num_qp_mother = np.array(num_qp_agg_species['mother'])
num_non_mother = np.array(num_non_agg_species['mother'])
qp_props_mother = num_qp_mother/(num_non_mother + num_qp_mother).astype('float') # one value for each tp

'''
num_qp = np.array(num_qp_agg_species['mother'] + num_qp_agg_species['infant'])
num_non = np.array(num_non_agg_species['mother'] + num_non_agg_species['infant'])
qp_props = num_qp/(num_non + num_qp).astype('float') # one value for each tp
'''

xticks = np.arange(len(labels))

fig, ax = plt.subplots(3, 1, figsize=(18, 12), sharex=True)
fig.subplots_adjust(wspace=0.02, hspace=0.02)

ax[0].boxplot(alpha_divs)
ax[0].set_ylabel("Shannon alpha diversity\nper sample")
ax[0].set_title("Alpha diversity, proportion of non-QP samples, polymorphism by timepoint (infants exclude Olm)")
ax[0].axvline(5.5, color='gray', linestyle='-')

# QP ===================================================================

ax[1].plot(xticks[:5] + 1, 1-qp_props_mother, 'r.-')
ax[1].plot(xticks[5:] + 1, 1-qp_props_infant, 'g.-')
ax[1].set_ylabel("Prop. of high coverage samples\nwhich are non-QP")
ax[1].axvline(5.5, color='gray', linestyle='-')
'''
ax[1].boxplot(prop_nonqp_sample['mother'], positions=(xticks[:5] + 1))
ax[1].boxplot(prop_nonqp_sample['infant'], positions=(xticks[5:] + 1))
ax[1].set_ylabel("Prop. of high coverage species in a sample\nwhich are non-QP")
ax[1].axvline(5.5, color='gray', linestyle='-')
'''
# ======================================================================

bplot_s = ax[2].boxplot(seg_pSs, patch_artist=True, widths=0.16, positions=(xticks+0.75))
for patch in bplot_s['boxes']:
	patch.set_facecolor('pink')

bplot_c = ax[2].boxplot(common_pSs, patch_artist=True, widths=0.16, positions=(xticks+1))
for patch in bplot_c['boxes']:
	patch.set_facecolor('lightblue')

bplot_r = ax[2].boxplot(rare_pSs, patch_artist=True, widths=0.16, positions=(xticks+1.25))
for patch in bplot_r['boxes']:
	patch.set_facecolor('lightgreen')

ax[2].set_yscale('log')
ax[2].set_ylabel("Polymorphism (pS)\nper sample-species pair")
ax[2].set_xlim((0, 24))
ax[2].set_xticks(xticks + 1)
ax[2].set_xticklabels(labels, fontsize=11)
ax[2].set_xlabel("Timepoint")
ax[2].axvline(5.5, color='gray', linestyle='-')

ax[1].legend([bp["boxes"][0] for bp in [bplot_s, bplot_c, bplot_r]], ['any (MAF >0)', 'common (MAF >0.2)', 'rare (MAF <0.05)'], loc='lower right')

fig.savefig("%s/alpha_div-prop_qp-polymorphism_over_time_v3.pdf" % config.analysis_directory, bbox_inches='tight')
