import matplotlib
matplotlib.use('Agg')
from utils import sample_utils as su, config, parse_midas_data
import sys, numpy
from math import log10,ceil
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# =================================================================
# Load pickled pN/pS data (from pickle_pnps.py)

type = sys.argv[1] # common, rare, etc.
if type == 'common':
	type_desc = '0.2-0.5'
elif type == 'rare':
	type_desc = '0-0.05 (excl. 0)'

import pickle
adir = config.analysis_directory
pdir = "%s/pickles" % config.data_directory

within_total_sites_nonQP = pickle.load(open("%s/within_total_sites_%s_nonQP.pkl" % (pdir, type), 'rb'))
within_total_sites_QP = pickle.load(open("%s/within_total_sites_%s_QP.pkl" % (pdir, type), 'rb'))
pNpS_nonQP = pickle.load(open("%s/pNpS_%s_nonQP.pkl" % (pdir, type), 'rb'))
pNpS_QP = pickle.load(open("%s/pNpS_%s_QP.pkl" % (pdir, type), 'rb'))
# =================================================================

# Timepoint-samples dictionary
# Exclude Olm (preterm infants)
tp_sample_dict, ordered_tps = su.get_mi_tp_sample_dict(exclude_cohorts = ['olm'], binned = True)
num_tp = len(tp_sample_dict)

# Colormap
colormap = cm.get_cmap('cool', num_tp)
colors = [colormap(x) for x in numpy.array([x for x in range(0,num_tp)]) / float(num_tp)]

# =================================================================
# Scatterplot of pN/pS vs. pS values per sample-species pair, QP
# Color-coded by timepoint
# =================================================================

fig, ax = plt.subplots(figsize=(10,8))
ax.set_xscale('log')
ax.set_ylim((0, 3.5)) # Missing a weird pN/pS = 7 point...

color_i = 0

for tp in ordered_tps:
	all_pNpS, all_pS = [], []
	for sample in tp_sample_dict[tp]:
		for species in within_total_sites_QP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_QP[sample][species]
			# Fraction of all nonsynonymous sites with minor allele frequency > 0.05
			pN = (w1D*1.0 + 1.0)/(t1D + 1.0)		
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			pS = (w4D*1.0 + 1.0)/(t4D + 1.0)
			all_pS.append(pS)
			all_pNpS.append(pN/pS)
	
	ax.plot(all_pS, all_pNpS, '.', color = colors[color_i], label=tp, markersize=3)	
	color_i += 1

ax.set_title("%s pN/pS per sample-species pair (QP samples only)" % type)
ax.set_xlabel("pS: Fraction of all synonymous sites with minor allele frequency %s" % type_desc)
ax.set_ylabel("pN/pS")
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

fig.tight_layout()
fig.savefig("%s/pNpS_per_sample_species_pair_%s_QP.pdf" % (adir, type), bbox_inches='tight')

# =================================================================
# Scatterplot of pN/pS vs. pS values per sample-species pair, nonQP
# Color-coded by timepoint
# =================================================================

fig, ax = plt.subplots(figsize=(8,8))
ax.set_xscale('log')
ax.set_ylim((0, 3.5)) # Missing a weird pN/pS = 7 point again...

color_i = 0

for tp in ordered_tps:
	all_pNpS, all_pS = [], []
	for sample in tp_sample_dict[tp]:
		for species in within_total_sites_nonQP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_nonQP[sample][species]
			# Fraction of all nonsynonymous sites with minor allele frequency > 0.05
			pN = (w1D*1.0 + 1.0)/(t1D + 1.0)		
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			pS = (w4D*1.0 + 1.0)/(t4D + 1.0)
			all_pS.append(pS)
			all_pNpS.append(pN/pS)
	
	ax.plot(all_pS, all_pNpS, '.', color = colors[color_i], label=tp)	
	color_i += 1

ax.set_title("%s pN/pS per sample-species pair (non-QP samples only)" % type)
ax.set_xlabel("pS: Fraction of all synonymous sites with minor allele frequency %s" % type_desc)
ax.set_ylabel("pN/pS")
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')

fig.savefig("%s/pNpS_per_sample_species_pair_%s_nonQP.pdf" % (adir, type), bbox_inches='tight')

# =================================================================
# Boxplots for pS per sample-species pair (no aggregation)
# =================================================================

fig, ax = plt.subplots(1, 2, figsize=(25,12), sharey=True)

ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_ylabel("pS: Fraction of all synonymous sites\nwith minor allele frequency %s" % type_desc)

# QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pS = []
	for sample in tp_sample_dict[tp]:
		for species in within_total_sites_QP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_QP[sample][species]
			# Fraction of all nonsynonymous sites with minor allele frequency > 0.05
			pN = (w1D*1.0 + 1.0)/(t1D + 1.0)		
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			pS = (w4D*1.0 + 1.0)/(t4D + 1.0)
			all_pS.append(pS)
	label = "%s\n(n=%i)" % (tp, len(all_pS))
	all_tp_labels.append(label)
	all_tp_pS.append(all_pS)

ax[0].boxplot(all_tp_pS)
ax[0].set_xticklabels(all_tp_labels)
ax[0].set_title("%s pS per sample-species pair (QP samples only)" % type)

# non-QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pS = []
	for sample in tp_sample_dict[tp]:
		for species in within_total_sites_nonQP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_nonQP[sample][species]
			# Fraction of all nonsynonymous sites with minor allele frequency > 0.05
			pN = (w1D*1.0 + 1.0)/(t1D + 1.0)		
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			pS = (w4D*1.0 + 1.0)/(t4D + 1.0)
			all_pS.append(pS)
	label = "%s\n(n=%i)" % (tp, len(all_pS))
	all_tp_labels.append(label)
	all_tp_pS.append(all_pS)

ax[1].boxplot(all_tp_pS)
ax[1].set_xticklabels(all_tp_labels)
ax[1].set_title("%s pS per sample-species pair (non-QP samples only)" % type)

fig.tight_layout()
fig.savefig("%s/pS_per_sample_species_pair_boxplots_%s.png" % (adir, type), bbox_inches="tight")

# =================================================================
# Boxplots for pS per sample (averaged over species)
# Color-coded by timepoint
# =================================================================

fig, ax = plt.subplots(1, 2, figsize=(25,12), sharey=True)

ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_ylabel("pS: Fraction of all synonymous sites\nwith minor allele frequency %s" % type_desc)

# QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pS = []
	for sample in tp_sample_dict[tp]:
		sample_w4D = 0
		sample_t4D = 0
		for species in within_total_sites_QP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_QP[sample][species]
			sample_w4D += w4D
			sample_t4D += t4D
		if sample_t4D != 0:
			pS = (sample_w4D*1.0 + 1.0)/(sample_t4D + 1.0)
			all_pS.append(pS)
	label = "%s\n(n=%i)" % (tp, len(all_pS))
	all_tp_labels.append(label)
	all_tp_pS.append(all_pS)

ax[0].boxplot(all_tp_pS)
ax[0].set_xticklabels(all_tp_labels)
ax[0].set_title("%s pS per sample (QP samples only)" % type)

# non-QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pS = []
	for sample in tp_sample_dict[tp]:
		sample_w4D = 0
		sample_t4D = 0
		for species in within_total_sites_nonQP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_nonQP[sample][species]
			sample_w4D += w4D
			sample_t4D += t4D
		if sample_t4D != 0:
			pS = (sample_w4D*1.0 + 1.0)/(sample_t4D + 1.0)
			all_pS.append(pS)
	label = "%s\n(n=%i)" % (tp, len(all_pS))
	all_tp_labels.append(label)
	all_tp_pS.append(all_pS)

ax[1].boxplot(all_tp_pS)
ax[1].set_xticklabels(all_tp_labels)
ax[1].set_title("%s pS per sample (non-QP samples only)" % type)

fig.tight_layout()
fig.savefig("%s/pS_per_sample_boxplots_%s.png" % (adir, type), bbox_inches="tight")

# =================================================================
# Boxplots for pN/pS per sample-species pair (no aggregation)
# Color-coded by timepoint
# =================================================================

fig, ax = plt.subplots(1, 2, figsize=(25,18), sharey=True)

# ax[0].set_yscale('log')
# ax[1].set_yscale('log')
# ax1.set_xticks([20, 200, 500])
# ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# ax[0].set_ylim((0.0001, 1))
if type == 'rare':
	ax[0].set_ylim((0, 3.5)) # Missing a weird pN/pS = 7 point...
ax[0].set_ylabel("pN/pS")

# QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pNpS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pNpS = []
	for sample in tp_sample_dict[tp]:
		for species in within_total_sites_QP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_QP[sample][species]
			# Fraction of all nonsynonymous sites with minor allele frequency > 0.05
			pN = (w1D*1.0 + 1.0)/(t1D + 1.0)		
			# Fraction of all synonymous sites with minor allele frequency > 0.05
			pS = (w4D*1.0 + 1.0)/(t4D + 1.0)
			all_pNpS.append(pN/pS)
	label = "%s\n(n=%i)" % (tp, len(all_pNpS))
	all_tp_labels.append(label)
	all_tp_pNpS.append(all_pNpS)

ax[0].boxplot(all_tp_pNpS)
ax[0].set_xticklabels(all_tp_labels)
ax[0].set_title("%s pN/pS per sample-species pair (QP samples only)" % type)

# non-QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pNpS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pNpS = []
	for sample in tp_sample_dict[tp]:
		for species in within_total_sites_nonQP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_nonQP[sample][species]
			# Fraction of all nonsynonymous sites with minor allele frequency in specified range
			pN = (w1D*1.0 + 1.0)/(t1D + 1.0)
			# Fraction of all synonymous sites with minor allele frequency in specified range
			pS = (w4D*1.0 + 1.0)/(t4D + 1.0)
			all_pNpS.append(pN/pS)
	label = "%s\n(n=%i)" % (tp, len(all_pNpS))
	all_tp_labels.append(label)
	all_tp_pNpS.append(all_pNpS)

ax[1].boxplot(all_tp_pNpS)
ax[1].set_xticklabels(all_tp_labels)
ax[1].set_title("%s pN/pS per sample-species pair (non-QP samples only)" % type)

fig.tight_layout()
fig.savefig("%s/pNpS_per_sample_species_pair_boxplots_%s.png" % (adir, type), bbox_inches="tight")

# =================================================================
# Boxplots for pN/pS per sample (averaged over species)
# Color-coded by timepoint
# =================================================================

fig, ax = plt.subplots(1, 2, figsize=(25,18), sharey=True)

# ax[0].set_yscale('log')
# ax[1].set_yscale('log')
if type == 'rare':
	ax[0].set_ylim((0, 3.5)) # Missing a weird pN/pS = 7 point...
# ax[0].set_ylim((0.0001, 1))
ax[0].set_ylabel("pN/pS")

# QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pNpS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pNpS = []
	for sample in tp_sample_dict[tp]:
		sample_w1D, sample_t1D, sample_w4D, sample_t4D = 0, 0, 0, 0
		for species in within_total_sites_QP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_QP[sample][species]
			sample_w1D += w1D
			sample_t1D += t1D
			sample_w4D += w4D
			sample_t4D += t4D
		if sample_t4D != 0 and sample_t1D != 0:
			pN = (sample_w1D*1.0 + 1.0)/(sample_t1D + 1.0)
			pS = (sample_w4D*1.0 + 1.0)/(sample_t4D + 1.0)
			all_pNpS.append(pN/pS)
	label = "%s\n(n=%i)" % (tp, len(all_pNpS))
	all_tp_labels.append(label)
	all_tp_pNpS.append(all_pNpS)

ax[0].boxplot(all_tp_pNpS)
ax[0].set_xticklabels(all_tp_labels)
ax[0].set_title("%s pN/pS per sample (QP samples only)" % type)

# non-QP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_tp_pNpS, all_tp_labels = [], []

for tp in ordered_tps:
	all_pNpS = []
	for sample in tp_sample_dict[tp]:
		sample_w1D, sample_t1D, sample_w4D, sample_t4D = 0, 0, 0, 0
		for species in within_total_sites_nonQP[sample]:
			w1D, t1D, w4D, t4D = within_total_sites_nonQP[sample][species]
			sample_w1D += w1D
			sample_t1D += t1D
			sample_w4D += w4D
			sample_t4D += t4D
		if sample_t4D != 0 and sample_t1D != 0:
			pN = (sample_w1D*1.0 + 1.0)/(sample_t1D + 1.0)
			pS = (sample_w4D*1.0 + 1.0)/(sample_t4D + 1.0)
			all_pNpS.append(pN/pS)
	label = "%s\n(n=%i)" % (tp, len(all_pNpS))
	all_tp_labels.append(label)
	all_tp_pNpS.append(all_pNpS)

ax[1].boxplot(all_tp_pNpS)
ax[1].set_xticklabels(all_tp_labels)
ax[1].set_title("%s pN/pS per sample (non-QP samples only)" % type)

fig.tight_layout()
fig.savefig("%s/pNpS_per_sample_boxplots_%s.png" % (adir, type), bbox_inches="tight")
