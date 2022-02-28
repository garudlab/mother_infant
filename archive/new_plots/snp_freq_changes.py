from utils import temporal_changes_utils, config, sample_utils as su, parse_midas_data
from collections import defaultdict
import numpy as np
import pickle
import sys

# Minimum coverage to consider sites
min_coverage = int(sys.argv[1]) # 20, 50, 75, 100

# Look at all allele frequency changes from very lax 0.45 -> 0.55
# I.e. frequency change of >= 10%
lower_threshold = 0.45
upper_threshold = 0.55

# Use alternate frequency and minimum coverage thresholds
# for determining whether modification occurred
min_coverage_mod = 20
lower_threshold_mod = 0.2
upper_threshold_mod = 0.8

cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'olm', 'hmp']

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.parse_sample_subject_map()
all_samples = sample_country_map.keys()
sys.stderr.write("Done!\n")

# Samples by cohort
mother_samples = su.get_sample_names('mother','all')
infant_samples = su.get_sample_names('infant','all')
hmp_samples = su.get_sample_names('hmp', 'all')

# Species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

if True:
	# Aggregate by species-cohort-timepoint pair
	# species -> cohort -> tp pair -> list of freq deltas
	snp_freq_deltas_dict = {s: {c: defaultdict(list) for c in cohorts} for s in good_species_list}
	
	# Aggregate by species-cohort-timepoint pair-subject
	# species -> cohort -> subject -> tp pair -> list of freq deltas
	snp_freq_deltas_subject_dict = {s: {c: defaultdict(dict) for c in cohorts} for s in good_species_list}
	
	# Map subject-timepoint pair to sample
	snp_subject_samples_dict = defaultdict(dict)

	for species_name in good_species_list:
			
			sys.stderr.write("\nProcessing %s...\n" % species_name)
			
			qp_sample_lists = {}
			for cohort in cohorts:
				qp_sample_lists[cohort] = sorted(su.load_qp_samples(su.get_sample_names(cohort,'all'), species_name)['qp'])
			
			sys.stderr.write("Loading pre-computed temporal changes...\n")
			temporal_change_map = temporal_changes_utils.load_temporal_change_map(species_name, min_coverage)
			temporal_change_map_mod = temporal_changes_utils.load_temporal_change_map(species_name, min_coverage_mod)
			sys.stderr.write("Done!\n")
			
			### Now loop over different cohorts
			for cohort in cohorts:		
				desired_samples = qp_sample_lists[cohort]
				
				same_subject_idxs = su.calculate_mi_ordered_same_subject_pairs(sample_order_map, desired_samples, within_host_type='consecutive', one_per_mi_pair=False)
				
				for i, j in zip(same_subject_idxs[0], same_subject_idxs[1]):
					
					sample_i = desired_samples[i] 
					sample_j = desired_samples[j]
					
					
					subject = sample_subject_map[sample_i]
					if subject.split('-')[0] != sample_subject_map[sample_j].split('-')[0]:
						print("Huh??? Subject doesn't match")
						continue
					
					tp_pair = su.sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples)
					
					snp_subject_samples_dict[subject][tp_pair] = (sample_i, sample_j)
					
					_, _, mutations_full, reversions_full = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map_mod, sample_i, sample_j, lower_threshold=lower_threshold_mod, upper_threshold=upper_threshold_mod)
					
					all_snp_changes_full = mutations_full + reversions_full
					
					L, perr, mutations, reversions = temporal_changes_utils.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold=lower_threshold, upper_threshold=upper_threshold)
					
					all_snp_changes = mutations + reversions
					
					# Skip replacements and things with high error rates
					if perr >= 0.5 or len(all_snp_changes_full) > 100:
						continue
					
					snp_freq_deltas_subject_dict[species_name][cohort][subject][tp_pair] = []
					
					for snp_change in all_snp_changes:
						gene_id, contig, position, variant_type, A1, D1, A2, D2 = snp_change
						freq_delta = abs((float(A1)/D1) - (float(A2)/D2))
						snp_freq_deltas_dict[species_name][cohort][tp_pair].append(freq_delta)
						snp_freq_deltas_subject_dict[species_name][cohort][subject][tp_pair].append(freq_delta)
	
	pickle.dump(snp_subject_samples_dict, open("%s/pickles/snp_subject_samples_dict.pkl" % (config.data_directory), "wb"))
	pickle.dump(snp_freq_deltas_dict, open("%s/pickles/snp_freq_deltas_alt/snp_freq_deltas_cov_%i_dict.pkl" % (config.data_directory, min_coverage), "wb"))
	pickle.dump(snp_freq_deltas_subject_dict, open("%s/pickles/snp_freq_deltas_alt/snp_freq_deltas_subject_cov_%i_dict.pkl" % (config.data_directory, min_coverage), "wb"))

snp_freq_deltas_dict = pickle.load(open("%s/pickles/snp_freq_deltas_alt/snp_freq_deltas_cov_%i_dict.pkl" % (config.data_directory, min_coverage), 'r'))
snp_freq_deltas_subject_dict = pickle.load(open("%s/pickles/snp_freq_deltas_alt/snp_freq_deltas_subject_cov_%i_dict.pkl" % (config.data_directory, min_coverage), "rb"))

def tp_to_bin(tp_pair):
	tpa, tpb = tp_pair
	order_a = float(tpa[1:])
	order_b = float(tpb[1:])
	o1, o2 = (order_a, order_b) if (order_a <= order_b) else (order_b, order_a)
	
	# Custom bins
	tp_bins = {(0, 7.5): "Week 1", (7.5, 14.5): "Week 2", (14.5, 31): "Week 3-4", (31, 61): "month 2", (61, 122): "month 3-4", (122, 244): "month 5-8", (244, 999): "month 9+"}
	
	for start_inc, end_exc in tp_bins:
		if o1 >= start_inc and o1 < end_exc:
			bin1 = tp_bins[(start_inc, end_exc)]
		if o2 >= start_inc and o2 < end_exc:
			bin2 = tp_bins[(start_inc, end_exc)]
	
	tp_cat = 'MI' if (tpa[0], tpb[0]) == ('I', 'M') else tpa[0]+tpb[0]
	string = tp_cat + ": " + bin1 + " > " + bin2
	return string

# Aggregate data by cohort, tp pair, and species

snp_freq_deltas_by_cohort = defaultdict(list)
snp_freq_deltas_by_tp_pair = defaultdict(list)
snp_freq_deltas_by_tp_type = defaultdict(list)
snp_freq_deltas_by_species = defaultdict(list)
snp_freq_deltas_by_tp_bin = defaultdict(list)

for species in snp_freq_deltas_dict:
	for cohort in snp_freq_deltas_dict[species]:
		for tp_pair in snp_freq_deltas_dict[species][cohort]:
			freq_deltas = snp_freq_deltas_dict[species][cohort][tp_pair]
			snp_freq_deltas_by_cohort[cohort] += freq_deltas
			snp_freq_deltas_by_tp_pair[tp_pair] += freq_deltas
			snp_freq_deltas_by_species[species] += freq_deltas
			
			tpa, tpb = tp_pair
			tp_type = tpa[0] + tpb[0]
			if tp_type == 'IM':
				tp_type = 'MI'
			snp_freq_deltas_by_tp_type[tp_type] += freq_deltas
			
			tp_bin = tp_to_bin(tp_pair)
			snp_freq_deltas_by_tp_bin[tp_bin] += freq_deltas

# Get species with the most SNP changes (any cohort, tp pair)
num_snp_changes = []
all_species = []

for species in snp_freq_deltas_by_species:
	num_snp_changes.append(len(snp_freq_deltas_by_species[species]))
	all_species.append(species)

num_snp_changes = np.array(num_snp_changes)
all_species = np.array(all_species)

sorted_all_species, _ = zip(*sorted(zip(all_species, num_snp_changes), key=lambda x: x[1], reverse=True))

big_species = list(sorted_all_species)[:11]

# More aggregation...

snp_freq_deltas_by_species_tp_type = defaultdict(list)
snp_freq_deltas_by_cohort_tp_type = defaultdict(list)

for species in snp_freq_deltas_dict:
	for cohort in snp_freq_deltas_dict[species]:
		for tp_pair in snp_freq_deltas_dict[species][cohort]:
			tpa, tpb = tp_pair
			tp_type = tpa[0] + tpb[0]
			if tp_type == 'IM':
				tp_type = 'MI'
			freq_deltas = snp_freq_deltas_dict[species][cohort][tp_pair]
			if species in big_species:
				snp_freq_deltas_by_species_tp_type[(species, tp_type)] += freq_deltas
			snp_freq_deltas_by_cohort_tp_type[(cohort, tp_type)] += freq_deltas

# Plot time!

import os

plot_dir = "%s/snp_freq_changes/cov_%i/" % (config.analysis_directory, min_coverage)
os.system('mkdir -p %s' % plot_dir)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as stats

# Plot histogram of SNP freq deltas by cohort
fig_cohort, ax_cohort = plt.subplots(1, 6, figsize=(28,6))

for i in range(len(cohorts)):
	cohort = cohorts[i]
	freq_deltas = snp_freq_deltas_by_cohort[cohort]
	# ax_cohort[i].hist(freq_deltas, bins=len(set(freq_deltas)))
	ax_cohort[i].hist(freq_deltas, bins=25)
	ax_cohort[i].set_ylabel('Number of SNP changes')
	ax_cohort[i].set_xlabel('Frequency difference')
	ax_cohort[i].set_title('Freq deltas, %s, n=%i' % (cohort, len(freq_deltas)))

sys.stderr.write("Saving figures...\n")
fig_cohort.savefig('%s/freq_deltas_by_cohort.pdf' % (plot_dir),bbox_inches='tight')

# Plot histogram of SNP freq deltas by cohort, but overlaid
fig_cohort_o, ax_cohort_o = plt.subplots(1, 1)

for i in range(len(cohorts)):
	cohort = cohorts[i]
	freq_deltas = snp_freq_deltas_by_cohort[cohort]
	bins = 40
	weights = np.ones_like(freq_deltas)/float(len(freq_deltas))
	# ax_cohort_o.hist(freq_deltas, bins, weights = weights, alpha = 0.3, label = cohort + " (n=%i)" % len(freq_deltas))
	n,x,_ = plt.hist(freq_deltas, bins, weights=weights, histtype=u'step', alpha = 0.0)
	bin_centers = 0.5*(x[1:]+x[:-1])
	ax_cohort_o.plot(bin_centers, n, alpha = 0.5, label = cohort + " (n=%i)" % len(freq_deltas))

# ax_cohort_o.set_ylabel('Density')
ax_cohort_o.set_xlabel('Frequency difference')
ax_cohort_o.set_title('Histograms for freq deltas by cohort')
ax_cohort_o.legend()

sys.stderr.write("Saving figures...\n")
fig_cohort_o.savefig('%s/freq_deltas_by_cohort_overlay.png' % (plot_dir),bbox_inches='tight')

# Plot histogram of SNP freq deltas by species, but overlaid
fig_species_o, ax_species_o = plt.subplots(1, 1)

for i in range(len(good_species_list)):
	species = good_species_list[i]
	freq_deltas = snp_freq_deltas_by_species[species]
	if len(freq_deltas) < 80:
		continue
	bins = 10
	weights = np.ones_like(freq_deltas)/float(len(freq_deltas))
	n,x,_ = plt.hist(freq_deltas, bins, weights=weights, histtype=u'step', alpha = 0.0)
	bin_centers = 0.5*(x[1:]+x[:-1])
	ax_species_o.plot(bin_centers, n, alpha = 0.5, label = species + " (n=%i)" % len(freq_deltas))

# ax_cohort_o.set_ylabel('Density')
ax_species_o.set_xlabel('Frequency difference')
ax_species_o.set_title('Histograms for freq deltas by species')

box = ax_species_o.get_position()
ax_species_o.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax_species_o.legend(loc='center left', bbox_to_anchor=(1, 0.5))

sys.stderr.write("Saving figures...\n")
fig_species_o.savefig('%s/freq_deltas_by_species_overlay.png' % (plot_dir),bbox_inches='tight')

# Plot histogram of SNP freq deltas by tp type, but overlaid
fig_tp_o, ax_tp_o = plt.subplots(1, 1)

for i in range(len(snp_freq_deltas_by_tp_type.keys())):
	tp_type = snp_freq_deltas_by_tp_type.keys()[i]	
	freq_deltas = snp_freq_deltas_by_tp_type[tp_type]
	bins = 20
	if len(freq_deltas) <= 1:
		continue
	
	weights = np.ones_like(freq_deltas)/float(len(freq_deltas))
	n,x,_ = plt.hist(freq_deltas, bins, weights=weights, histtype=u'step', alpha = 0.0)
	bin_centers = 0.5*(x[1:]+x[:-1])
	ax_tp_o.plot(bin_centers, n, alpha = 0.5, label = tp_type + " (n=%i)" % len(freq_deltas))

# ax_cohort_o.set_ylabel('Density')
ax_tp_o.set_xlabel('Frequency difference')
ax_tp_o.set_title('Histograms for freq deltas by tp type')
ax_tp_o.legend()

sys.stderr.write("Saving figures...\n")
fig_tp_o.savefig('%s/freq_deltas_by_tp_type_overlay.png' % (plot_dir),bbox_inches='tight')

# ====================================================================
# Plot histogram of SNP freq deltas by tp type, stacked
# ====================================================================

tp_type_to_string = {'II': 'Infant-Infant', 'MI': 'Mother-Infant', 'MM': 'Mother-Mother', 'AA': 'Adult-Adult'}

fig_tp, ax_tp = plt.subplots(4, 1, figsize=(6,8), sharex=True)

ordered_tp_type = ['II', 'MI', 'MM', 'AA']

for i in range(len(ordered_tp_type)):
	tp_type = ordered_tp_type[i]
	freq_deltas = snp_freq_deltas_by_tp_type[tp_type]
	bins = 20
	
	if len(freq_deltas) <= 1:
		continue
	
	n,x,_ = plt.hist(freq_deltas, bins, histtype=u'step', alpha = 0.0)
	bin_centers = 0.5*(x[1:]+x[:-1])
	ax_tp[i].plot(bin_centers, n, 'b-')
	if i == 0:
		ax_tp[i].set_title('Histograms for freq deltas by tp type\n' +
		tp_type + " (n=%i)" % len(freq_deltas))
	else:
		ax_tp[i].set_title(tp_type + " (n=%i)" % len(freq_deltas))

# ax_cohort_o.set_ylabel('Density')
ax_tp[3].set_xlabel('Alt allele frequency difference')
ax_tp[3].set_ylim((0, (int(max(n)/100)+1)*100))

sys.stderr.write("Saving figures...\n")
fig_tp.savefig('%s/freq_deltas_by_tp_type_stacked.png' % (plot_dir),bbox_inches='tight')

# ====================================================================
# Plot histogram of SNP freq deltas by tp bin, stacked
# ====================================================================

fig_tp, ax_tp = plt.subplots(25, 1, figsize=(6,40), sharex=True)

ordered_tp_bins = sorted(snp_freq_deltas_by_tp_bin.keys())

for i in range(len(ordered_tp_bins)):
	tp_bin = ordered_tp_bins[i]
	freq_deltas = snp_freq_deltas_by_tp_bin[tp_bin]
	bins = 20
	
	if len(freq_deltas) <= 1:
		continue
	
	n,x,_ = plt.hist(freq_deltas, bins, histtype=u'step', alpha = 0.0)
	bin_centers = 0.5*(x[1:]+x[:-1])
	ax_tp[i].plot(bin_centers, n, 'b-')
	if i == 0:
		ax_tp[i].set_title('Histograms for freq deltas by tp bin\n' +
		tp_bin + " (n=%i)" % len(freq_deltas))
	else:
		ax_tp[i].set_title(tp_bin + " (n=%i)" % len(freq_deltas))

# ax_cohort_o.set_ylabel('Density')
ax_tp[24].set_xlabel('Alt allele frequency difference')
ax_tp[24].set_ylim((0, (int(max(n)/10)+1)*10))

sys.stderr.write("Saving figures...\n")
fig_tp.savefig('%s/freq_deltas_by_tp_bin_stacked.png' % (plot_dir),bbox_inches='tight')

# ====================================================================
# Plot SNP freq deltas by tp type x cohort
# ====================================================================

tp_type_to_string = {'II': 'Infant-Infant', 'MI': 'Mother-Infant', 'MM': 'Mother-Mother', 'AA': 'Adult-Adult'}

cohorts_reformat = ['Backhed', 'Ferretti', 'Yassour', 'Shao', 'Olm', 'HMP']
tp_types = ['II', 'MI', 'MM', 'AA']

fig_tp, ax_tp = plt.subplots(4, 6, figsize=(20, 12), sharex=True)

for tp_i in range(len(tp_types)):
	for ch_i in range(len(cohorts)):
		
		tp_type, cohort = tp_types[tp_i], cohorts[ch_i]
		
		if tp_i == 0 and ch_i == 0:
			ax_tp[tp_i][ch_i].set_title(cohorts_reformat[ch_i] + '\n' + tp_type_to_string[tp_type])
		elif tp_i == 0:
			ax_tp[tp_i][ch_i].set_title(cohorts_reformat[ch_i], loc='center')
		elif ch_i == 0:
			ax_tp[tp_i][ch_i].set_title(tp_type_to_string[tp_type], loc='left')
		
		if (cohort, tp_type) not in snp_freq_deltas_by_cohort_tp_type:
			continue
		freq_deltas = snp_freq_deltas_by_cohort_tp_type[(cohort, tp_type)]
		if len(freq_deltas) <= 1:
			continue
		
		bins=20
		n,x,_ = ax_tp[tp_i][ch_i].hist(freq_deltas, bins, histtype=u'step', alpha = 0.0)
		bin_centers = 0.5*(x[1:]+x[:-1])
		ax_tp[tp_i][ch_i].plot(bin_centers, n, 'b-', label=("(n=%i)" % len(freq_deltas)))
		ax_tp[tp_i][ch_i].legend(fontsize=8)

# ax_cohort_o.set_ylabel('Density')
ax_tp[3][1].set_xlabel('Alt allele frequency difference')

sys.stderr.write("Saving figures...\n")
fig_tp.savefig('%s/freq_deltas_by_cohort_and_tp_type.png' % (plot_dir),bbox_inches='tight')

# ====================================================================
# Plot SNP freq deltas by tp type x species
# ====================================================================

tp_type_to_string = {'II': 'Infant-Infant', 'MI': 'Mother-Infant', 'MM': 'Mother-Mother', 'AA': 'Adult-Adult'}

big_species_reformat = ['\n'.join(sname.split('_')) for sname in big_species]

tp_types = ['II', 'MI', 'MM', 'AA']

fig_tp, ax_tp = plt.subplots(4, 11, figsize=(28,10), sharex=True)

for tp_i in range(len(tp_types)):
	for sp_i in range(len(big_species)):
		
		tp_type, species = tp_types[tp_i], big_species[sp_i]
		
		if tp_i == 0 and sp_i == 0:
			ax_tp[tp_i][sp_i].set_title(big_species_reformat[sp_i] + '\n' + tp_type_to_string[tp_type])
		elif tp_i == 0:
			ax_tp[tp_i][sp_i].set_title(big_species_reformat[sp_i], loc='center')
		elif sp_i == 0:
			ax_tp[tp_i][sp_i].set_title(tp_type_to_string[tp_type], loc='left')
		
		if (species, tp_type) not in snp_freq_deltas_by_species_tp_type:
			continue
		freq_deltas = snp_freq_deltas_by_species_tp_type[(species, tp_type)]
		
		if len(freq_deltas) <= 1:
			continue
		
		# bins = 20
		n,x,_ = plt.hist(freq_deltas, histtype=u'step', alpha = 0.0)
		bin_centers = 0.5*(x[1:]+x[:-1])
		ax_tp[tp_i][sp_i].plot(bin_centers, n, 'b-', label=("(n=%i)" % len(freq_deltas)))
		ax_tp[tp_i][sp_i].legend(fontsize=8)		

# ax_cohort_o.set_ylabel('Density')
ax_tp[3][5].set_xlabel('Alt allele frequency difference')
ax_tp[3][10].set_ylim((0, 30))

sys.stderr.write("Saving figures...\n")
fig_tp.savefig('%s/freq_deltas_by_species_and_tp_type.png' % (plot_dir),bbox_inches='tight')
