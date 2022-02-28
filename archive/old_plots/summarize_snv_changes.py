import pickle
import config
import numpy as np
from collections import defaultdict

# Parameters
cohorts = ['backhed', 'ferretti', 'yassour', 'hmp']
event_types = ["no change", "modification", "replacement"]
modification_threshold = 100
replacement_threshold = 400

# Pickle
pooled_snp_pickle_full_fn = '%s/pickles/pooled_snp_change_full.pkl' % (config.data_directory)
pooled_snp_change_distribution_full = pickle.load(open(pooled_snp_pickle_full_fn, 'rb'))

pooled_snp_pickle_all_fn = '%s/pickles/pooled_snp_change_partial.pkl' % (config.data_directory)
pooled_snp_change_distribution_all = pickle.load(open(pooled_snp_pickle_all_fn, 'rb'))

# Manually fix full/all discrepancies
pooled_snp_change_distribution_all['yassour'][frozenset(['M1', 'M3'])].pop(23)
pooled_snp_change_distribution_all['yassour'][frozenset(['M1', 'M2'])].pop(2)

# Build partial sweep only distribution and percentages
pooled_snp_change_distribution_partial = {cohort: {} for cohort in cohorts}
pooled_snp_change_distribution_partial_percents = {cohort: {} for cohort in cohorts}

for cohort in cohorts:
	for tp in pooled_snp_change_distribution_all[cohort]:
		all_list = np.array(pooled_snp_change_distribution_all[cohort][tp])
		full_list = np.array(pooled_snp_change_distribution_full[cohort][tp])
		partial_list = list(all_list - full_list)
		pooled_snp_change_distribution_partial[cohort][tp] = partial_list
		pooled_snp_change_distribution_partial_percents[cohort][tp] = list(partial_list/all_list.astype(float))

# Outputs
full_sweep_summary = {event: defaultdict(int) for event in event_types}
all_sweep_summary = {event: defaultdict(int) for event in event_types}

full_sweep_counts = {event: defaultdict(int) for event in event_types}
partial_sweep_counts = {event: defaultdict(int) for event in event_types}

average_partial_percents = {event: defaultdict(list) for event in ['modification', 'replacement']}

mod_or_replace_not_sure = []

event_types_alt = ["no change", "mod:maj-partial", "mod:maj-full", "replacement"]
all_sweep_summary_alt = {event: defaultdict(int) for event in event_types_alt}

def tp_to_category(tp_pair):
	tpa, tpb = tp_pair
	return tpa[0]+tpb[0]

for cohort in cohorts:
	full_subset = pooled_snp_change_distribution_full[cohort]
	partial_subset = pooled_snp_change_distribution_partial[cohort]
	partial_percent_subset = pooled_snp_change_distribution_partial_percents[cohort]
	all_subset = pooled_snp_change_distribution_all[cohort]
	for tp_pair in full_subset.keys():
		cat = tp_to_category(tp_pair)
		# Full
		for snp_change_count in full_subset[tp_pair]:
			if snp_change_count == 0:
				full_sweep_summary["no change"][cat] += 1
				full_sweep_counts["no change"][cat] += snp_change_count
			elif snp_change_count <= modification_threshold:
				full_sweep_summary["modification"][cat] += 1
				full_sweep_counts["modification"][cat] += snp_change_count
			elif snp_change_count >= replacement_threshold:
				full_sweep_summary["replacement"][cat] += 1
				full_sweep_counts["replacement"][cat] += snp_change_count
			else:
				mod_or_replace_not_sure.append(snp_change_count)
		# All
		for snp_change_count, full_snp_change_count, partial_snp_change_count, partial_percent in zip(all_subset[tp_pair], full_subset[tp_pair], partial_subset[tp_pair], partial_percent_subset[tp_pair]):
			if snp_change_count == 0:
				all_sweep_summary_alt["no change"][cat] += 1
				all_sweep_summary["no change"][cat] += 1
				partial_sweep_counts["no change"][cat] += partial_snp_change_count
			elif snp_change_count <= modification_threshold:
				all_sweep_summary["modification"][cat] += 1
				partial_sweep_counts["modification"][cat] += partial_snp_change_count
				average_partial_percents["modification"][cat].append(partial_percent)
				if partial_snp_change_count/float(snp_change_count) > 0.5:
					all_sweep_summary_alt["mod:maj-partial"][cat] += 1
				else:
					all_sweep_summary_alt["mod:maj-full"][cat] += 1
			elif snp_change_count >= replacement_threshold:
				all_sweep_summary_alt["replacement"][cat] += 1
				all_sweep_summary["replacement"][cat] += 1
				partial_sweep_counts["replacement"][cat] += partial_snp_change_count
				average_partial_percents["replacement"][cat].append(partial_percent)
			else:
				mod_or_replace_not_sure.append(snp_change_count)

for event in ['modification', 'replacement']:
	for tp in average_partial_percents[event]:
		average_partial_percents[event][tp] = np.mean(average_partial_percents[event][tp])

# For debugging purposes
cohort = 'yassour'
for x in pooled_snp_change_distribution_all[cohort]:
	all_len = len(pooled_snp_change_distribution_all[cohort][x])
	full_len = len(pooled_snp_change_distribution_full[cohort][x])
	print("all: " + str(all_len) + "| full: " + str(full_len))
	if all_len != full_len:
		print(x)

def normalize(arr):
	return np.array(arr)/float(sum(arr))
