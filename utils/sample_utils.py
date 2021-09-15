import numpy as np
import config
import parse_midas_data
from collections import defaultdict

# ===========================================================================
# This module (parse_metadata) contains the following utilities:
# 
# 	sample_to_tp, sample_pair_to_tp_pair
# 
# 	is_same_mi_subject, is_mi_pair
# 
# 	get_sample_names, get_ferretti_sample_names_by_site, calculate_qp_samples
# 	calculate_ordered_same_subject_pairs, calculate_ordered_diff_subject_pairs
# 	calculate_ordered_same_sample_pairs, calculate_ordered_subject_pairs
# 	calculate_mi_ordered_same_subject_pairs, calculate_sample_idx_map
# 	apply_sample_index_map_to_indices, calculate_samples_in_different_subjects
# 	calculate_subject_pairs, calculate_sample_subject_matrix
# 
# 	parse_sample_metadata_map, filter_sample_metadata_map, extract_sample_metadata_map
# 	parse_sample_subject_map, parse_sample_country_map, parse_sample_continent_map
# 	parse_subject_sample_map, parse_sample_order_map
# 
# 	flatten_samples, flatten_subjects, parse_isolate_metadata_map
# 	list_of_isolates_and_mixtures, parse_merged_sample_names
# 	calculate_unique_samples
# ===========================================================================

# Generic helper functions

flatten = lambda l: [item for sublist in l for item in sublist]

def subject_and_tp_to_sample(subject, tp, subject_sample_map, sample_order_map):
	for sample in subject_sample_map[subject].keys():
		subject, cur_order = sample_order_map[sample]
		if float(tp[1:]) == cur_order:
			return sample
	print("Sample not found")
	return None
	
# ===========================================================================
# FUNCTIONS FOR REFORMATTING METADATA
# ===========================================================================

# ===========================================================================
# sample_pair_to_tp: converts sample to generic timepoint descriptor
# A for adult, M for mother, I for infant
# ===========================================================================

def sample_to_tp(sample, sample_order_map, hmp_samples, mother_samples):
	order_str = str(sample_order_map[sample][1])
	if sample in hmp_samples:
		return 'A' + order_str
	else: # If not HMP, assume mother or infant
		return ('M' if sample in mother_samples else 'I') + order_str

# ===========================================================================
# sample_pair_to_tp_pair: converts sample pair to general timepoint pair form
# example: frozenset('I1', 'I3')
# A for adult, M for mother, I for infant
# ===========================================================================

def sample_pair_to_tp_pair(sample_i, sample_j, sample_order_map, hmp_samples, mother_samples):
	if sample_i in hmp_samples:
		tp_i = 'A' + str(sample_order_map[sample_i][1])
	else: # If not HMP, assume mother or infant
		tp_i = ('M' if sample_i in mother_samples else 'I') + str(sample_order_map[sample_i][1])
	if sample_j in hmp_samples:
		tp_j = 'A' + str(sample_order_map[sample_j][1])
	else: # If not HMP, assume mother or infant
		tp_j = ('M' if sample_j in mother_samples else 'I') + str(sample_order_map[sample_j][1])
	tp_pair = frozenset((tp_i, tp_j))
	return tp_pair

# ===========================================================================
# FUNCTIONS FOR DEALING WITH SAME SUBJECT MOTHER-INFANT
# ===========================================================================

# ===========================================================================
# get_same_mi_pair_dict: returns mapping of mother/child subject to
# child/mother subject (only valid for mi cohorts excluding olm)
# ===========================================================================

def get_same_mi_pair_dict(subject_sample_map):
	same_mi_pair_dict = {}
	for subject1 in subject_sample_map:
		for subject2 in subject_sample_map:
			if subject2 != subject1 and subject2[:-2] == subject1[:-2]:
				same_mi_pair_dict[subject1] = subject2
	return same_mi_pair_dict

# ===========================================================================
# is_same_mi_subject: returns true if two samples are from same mother-infant pair
# (Includes two infant or two mother samples at different timepoints)
# ===========================================================================

def is_same_mi_subject(sample1, sample2, sample_subject_map):
	
	# Note: this only works on mother/infant datasets
	# Observe the following subject name format conventions
	# Backhed: 64-I/-M
	# Ferretti: CA_C10006IS
	# Yassour: M0333-M
	# Shao: B02752-I
	# Olm: S2_017
	
	subject1 = sample_subject_map[sample1]
	subject2 = sample_subject_map[sample2]
	
	return (subject1[:-2] == subject2[:-2])

# ===========================================================================
# is_mi_pair: returns true if one sample is mother and other is infant
# ===========================================================================

def is_mi_pair(sample1, sample2, mother_samples, infant_samples):
	
	return ((sample1 in mother_samples and sample2 in infant_samples) or 
	(sample1 in infant_samples and sample2 in mother_samples))

# ===========================================================================
# FUNCTIONS FOR DEALING WITH SAMPLE LISTS
# ===========================================================================

# ===========================================================================
# get_sample_names: returns list of samples in specified cohort/timepoint
# 
# cohort: one of 'HMP', 'Backhed', 'Ferretti', 'Yassour', Shao', 'Olm'
# (case insensitive),	'mother', 'infant', 'adult', 'all', or 'all-dict'
# note that 'dict' is (cohort->tp->samples)
# 
# timepoint:
# 	HMP - 1, 2, 3
#		Backhed - B, 4M, 12M, M
#		Ferretti - M0, I1, I2, I3, I4, I5
#		Yassour - MGest, MBirth, M3, CBirth, C14, C1, C2, C3
# 	Shao - Mother, Neonatal (4-21 days), Infancy
# 	Olm - NIH1, NIH2, NIH3, NIH4, Sloan2
#		or 'all' for every sample in a cohort
# ===========================================================================

def get_sample_names(cohort, timepoint = 'all', remove_c=True):
	
	sample_dict = {'hmp': defaultdict(set), 'backhed': defaultdict(set), \
								 'ferretti': defaultdict(set), 'yassour': defaultdict(set), \
								 'shao': defaultdict(set), 'olm': defaultdict(set)}
	
	# Note: we need the sample lists because not all samples in the metadata
	# files should be returned. The final ID is based on metadata file, however
	# (so -c suffixes for combined samples will be omitted.)
	
	metadata_dir = config.metadata_directory
	samples_dir = "%s/final_sample_lists" % metadata_dir
	
	# HMP - expect 469 samples
	samples_fpath = "%s/HMP1-2_samples.txt" % samples_dir
	hmp_samples = [line.strip() for line in open(samples_fpath, 'r')]
	hmp_samples = parse_merged_sample_names(hmp_samples) if remove_c else hmp_samples # Remove c's
	
	metadata_fpath = "%s/HMP1-2_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]
	
	for sample in metadata[1:]:
		_, sample_id, _, _, _, tp = sample
		if sample_id in hmp_samples:
			sample_dict['hmp'][tp].add(sample_id)
		if (sample_id+'c') in hmp_samples:
			sample_dict['hmp'][tp].add((sample_id+'c'))
	
	# Backhed - expect 391 samples
	samples_fpath = "%s/Backhed_samples.txt" % samples_dir
	backhed_samples = [line.strip() for line in open(samples_fpath, 'r')]
	
	metadata_fpath = "%s/Backhed_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]
	
	for sample in metadata[1:]:
		sample_id, _, tp = sample
		if sample_id in backhed_samples:
			sample_dict['backhed'][tp].add(sample_id)
	
	# Ferretti - expect 119 samples
	samples_fpath = "%s/Ferretti_samples.txt" % samples_dir
	ferretti_samples = [line.strip() for line in open(samples_fpath, 'r')]
	
	metadata_fpath = "%s/Ferretti_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]
	
	for sample in metadata[1:]:
		sample_id = sample[4]
		tp = sample[6][9] + sample[6][19]
		# Restrict to fecal samples
		body_site = sample[6][15:17]
		if body_site == 'FE' and sample_id in ferretti_samples:
			sample_dict['ferretti'][tp].add(sample_id)
	
	# Yassour - expect 286 samples
	samples_fpath = "%s/Yassour_samples.txt" % samples_dir
	yassour_samples = [line.strip() for line in open(samples_fpath, 'r')]
	
	metadata_fpath = "%s/Yassour_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]
	
	for sample in metadata[1:]:
		sample_id = sample[4]
		tp_raw = sample[7][6:].split(':')
		tp = tp_raw[0][0] + tp_raw[1].split()[0]
		if sample_id in yassour_samples:
			sample_dict['yassour'][tp].add(sample_id)
	
	# Shao - expect 1679 samples
	samples_fpath = "%s/Shao_samples.txt" % samples_dir
	shao_samples = [line.strip() for line in open(samples_fpath, 'r')]
	
	metadata_fpath = "%s/Shao_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]
	
	for sample in metadata[1:]:
		sample_id, _, _, neonatal_tp, tp_cat = sample
		if tp_cat == 'NA':
			continue # Ignore if order info not known
		tp_cat = 'Infancy' if neonatal_tp == 'Infancy' else tp_cat
		if sample_id in shao_samples:
			sample_dict['shao'][tp_cat].add(sample_id)
	
	# Olm - expect 898 samples
	olm_samples = []
	for campaign in ['NIH1', 'NIH2', 'NIH3', 'NIH4', 'Sloan2']:
		samples_fpath = "%s/Olm_%s_samples.txt" % (samples_dir, campaign)
		olm_sub_samples = [line.strip() for line in open(samples_fpath, 'r')]
		if remove_c:
			olm_samples += list(parse_merged_sample_names(olm_sub_samples)) # Remove c's
		else:
			olm_samples += list(olm_sub_samples)
	
	metadata_fpath = "%s/Olm_metadata.txt" % metadata_dir
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file]
	
	for sample in metadata[1:]:
		try:
			sample_id = sample[9]
			campaign = sample[3]
			if sample_id in olm_samples:
				sample_dict['olm'][campaign].add(sample_id)
			if (sample_id+'c') in olm_samples:
				sample_dict['olm'][campaign].add((sample_id+'c'))
		except:
			continue
	
	# Aggregate + convert sets to lists
	all_samples = []
	for lcohort in sample_dict:
		for ltp in sample_dict[lcohort]:
			sample_dict[lcohort][ltp] = list(sample_dict[lcohort][ltp])
			all_samples += sample_dict[lcohort][ltp]
	
	# All HMP samples
	hmp_samples = flatten([sample_dict['hmp'][htp] for htp in sample_dict['hmp']])
	
	# All mother and infant samples
	mother_samples = []
	mother_samples += sample_dict['backhed']['M']
	mother_samples += sample_dict['ferretti']['M0']
	mother_samples += sample_dict['shao']['Mother']
	mother_samples += flatten([sample_dict['yassour'][ytp] for ytp in ['MGest', 'MBirth', 'M3']])
	infant_samples = [s for s in all_samples if (s not in mother_samples and s not in hmp_samples)]
	
	general_cohort_dict = {'all': all_samples, 'mother': mother_samples, 'infant': infant_samples, 'adult': mother_samples + hmp_samples}
	
	cohort = cohort.lower() # for case insensitivity
	
	if cohort in general_cohort_dict: # all, mother, infant, adult
		return general_cohort_dict[cohort]
	elif cohort == 'all-dict': # all-dict
		return sample_dict
	elif timepoint == 'all': # specific cohort, all
		all_cohort_samples = []
		for ltp in sample_dict[cohort]:
			all_cohort_samples += sample_dict[cohort][ltp]
		return all_cohort_samples
	else: # specific cohort, specific timepoint
		return sample_dict[cohort][timepoint]

# ===========================================================================
# get_mi_tp_sample_dict: returns dictionary that maps each each timepoint to
# list of samples across ALL mother-infant cohorts
# 
# The following conventions are used to unify timepoints across cohorts:
# 	- Assume each "month" is 30.5 days long, so month 1 = day 30.5
#		- Assume each "week" is 7 days long, so week 2 = day 14
# 	- Time values are expressed in days relative to birth/delivery
# 	- Time values are rounded to the nearest day
# 
# If binned = True, use each day as timepoint. Otherwise, use custom bins
# of timepoints (days 0-6, weeks 1-3, months 1-11, year 1 for infants)
# 
# ===========================================================================

def get_mi_tp_sample_dict(exclude_cohorts = [], binned = False):
	
	sample_cohort_dict = get_sample_names('all-dict')
	sample_order_map = parse_sample_order_map()
	tp_sample_dict = {'mother': defaultdict(list), 'infant': defaultdict(list)}
	
	# Backhed timepoints:
	# M: Mother, delivery / 0-5 days after birth (median 2)
	# B: Infant, birth / 2-5 days after birth (median 3)
	# 12M: Infant, 4 months / 119-125 days after birth (median 122)
	# 4M: Infant, 12 months / 363-372 days after birth (median 366)
	
	if 'backhed' not in exclude_cohorts:
		backhed_tp_map = {'M': 2, 'B': 3, '4M': 122, '12M': 366}
		for btp in backhed_tp_map:
			stp = backhed_tp_map[btp]
			mother_or_infant = 'mother' if (btp == 'M') else 'infant'
			tp_sample_dict[mother_or_infant][stp] += sample_cohort_dict['backhed'][btp]
	
	# Ferretti timepoints:
	# M0: Mother, delivery
	# I1: Infant, 1 day
	# I2: Infant, 3 days
	# I3: Infant, 1 week
	# I4: Infant, 1 month
	# I5: Infant, 4 months
	
	if 'ferretti' not in exclude_cohorts:
		ferretti_tp_map = {'M0': 0, 'I1': 1, 'I2': 3, 'I3': 7, 'I4': 30, 'I5': 122}
		for ftp in ferretti_tp_map:
			stp = ferretti_tp_map[ftp]
			mother_or_infant = 'mother' if (ftp == 'M0') else 'infant'
			tp_sample_dict[mother_or_infant][stp] += sample_cohort_dict['ferretti'][ftp]
	
	# Yassour timepoints:
	# MGest: Mother, gestational week 27
	# MBirth: Mother, delivery
	# M3: Mother, 3 months post-delivery
	# CBirth: Infant, birth / meconium
	# C14: Infant, 2 weeks
	# C1: Infant, 1 week
	# C2: Infant, 2 months
	# C3: Infant, 3 months
	
	if 'yassour' not in exclude_cohorts:
		yassour_tp_map = {'MGest': -92, 'MBirth': 0, 'M3': 92, 'CBirth': 0, 'C14': 14, 'C1': 30, 'C2': 61, 'C3': 92}
		for ytp in yassour_tp_map:
			stp = yassour_tp_map[ytp]
			mother_or_infant = 'mother' if (ytp[0] == 'M') else 'infant'
			tp_sample_dict[mother_or_infant][stp] += sample_cohort_dict['yassour'][ytp]
	
	# Shao timepoints:
	# Mother: one timepoint at delivery
	# Neonatal: 4-21 days (various)
	# Infancy: 4-15 months (various)
	
	if 'shao' not in exclude_cohorts:
		shao_samples = get_sample_names('Shao')
		for sample in shao_samples:
			subject, order = sample_order_map[sample]
			order = int(round(order))
			mother_or_infant = 'mother' if (subject[-1] == 'M') else 'infant'
			tp_sample_dict[mother_or_infant][order].append(sample)
	
	# Olm timepoints:
	# 5-86 days
	
	if 'olm' not in exclude_cohorts:
		olm_samples = get_sample_names('Olm')
		for sample in olm_samples:
			subject, order = sample_order_map[sample]
			order = int(round(order))
			tp_sample_dict['infant'][order].append(sample)
	
	if binned == False:
		return tp_sample_dict
	
	# Custom timepoint binning (label -> lower bound of bin, inclusive)
	custom_infant_bins = {'birth': 0, '1 day': 1, '2 day': 2, '3 day': 3, '4 day': 4, '5 day': 5, '6 day': 6, \
	'1 wk': 7, '2 wk': 14, '3 wk': 21, '1 mon': 30, '2 mon': 61, '3 mon': 91, '4 mon': 122, '5 mon': 152, '6 mon': 183, \
	'7 mon': 213, '8 mon': 244, '9 mon': 274, '10 mon': 305, '11 mon': 335, '1 yr': 365}
	custom_infant_bins_inverted = {v: k for k, v in custom_infant_bins.items()}
	custom_infant_bin_bounds = sorted(custom_infant_bins.values(), reverse=True)	
	
	# Like tp_sample_dict but uses new bins instead of days
	binned_infant_tp_sample_dict = defaultdict(list)
	
	for day in tp_sample_dict['infant']:
		# Find the right bin for this day
		for bound in custom_infant_bin_bounds:
			if day >= bound:
				bin_label = custom_infant_bins_inverted[bound]
				binned_infant_tp_sample_dict[bin_label] += tp_sample_dict['infant'][day]
				break
	
	custom_infant_bin_labels = []
	for bound in sorted(custom_infant_bins.values()):
		bin_label = custom_infant_bins_inverted[bound]
		if bin_label in binned_infant_tp_sample_dict:
			custom_infant_bin_labels.append(bin_label)
	
	return {'mother': tp_sample_dict['mother'], 'infant': binned_infant_tp_sample_dict}, custom_infant_bin_labels

def get_mi_sample_day_dict(exclude_cohorts = [], binned = False):
	
	mi_tp_sample_dict = get_mi_tp_sample_dict(exclude_cohorts, binned)
	mi_sample_day_dict = {}
	
	for cat in ['infant', 'mother']:
		for day in mi_tp_sample_dict[cat]:
			for sample in mi_tp_sample_dict[cat][day]:
				mi_sample_day_dict[sample] = day
	
	return mi_sample_day_dict

# ===========================================================================
# get_ferretti_sample_names_by_site: get dictionary of samples by body site
# Site names: FE (stool), SA (oral cavity), SK (skin), VA (vagina)
# ===========================================================================

def get_ferretti_sample_names_by_site(body_site_code):
	
	metadata_fpath = config.metadata_directory + "PRJNA352475.txt" # Ferretti
	with open(metadata_fpath, 'r') as metadata_file:
		metadata = [row.strip().split('\t') for row in metadata_file.readlines()]
	
	site_dict = defaultdict(list)
	
	for sample in metadata[1:]:
		sample_name = sample[4]
		body_site = sample[6][-8:-6]
		site_dict[body_site].append(sample_name)
	
	return site_dict[body_site_code]

# ===========================================================================
# Output information on the timepoint pairs present in a sample set
# ===========================================================================

# TODO

# ===========================================================================
# calculate_qp_samples: returns QP status dictionary given samples, species
# 'qp', 'non-qp' [high coverage], 'low-coverage'
# ===========================================================================

def calculate_qp_samples(all_samples, species_name, prev_cohort='all'):
	
	import diversity_utils
	
	# list of samples that meet coverage criteria for this species
	highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
	
	# list of samples that meet QP criteria for this species
	haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name, prev_cohort=prev_cohort))
	
	qp_statuses = ['qp','non-qp','low-coverage']
	qp_sample_sets = {status: set() for status in qp_statuses}
	
	for sample in all_samples:
		if sample in haploid_samples:	# QP (must be high coverage)
			qp_sample_sets['qp'].add(sample)
		elif sample in highcoverage_samples:	# Non-QP (high coverage)
			qp_sample_sets['non-qp'].add(sample)
		else:	# Low coverage
			qp_sample_sets['low-coverage'].add(sample)
	
	return qp_sample_sets

# ===========================================================================
# load_qp_samples: returns QP status dictionary given samples, species
# 'qp', 'non-qp' [high coverage], 'low-coverage' (uses pickle)
# ===========================================================================

def load_qp_samples(desired_samples, species_name, prev_cohort='all', force_repickle = False):
	
	import pickle, os.path
	pickle_fn = "%s/pickles/qp_samples/%s_qp_sample_dict.pkl" % (config.data_directory, species_name)
	
	if force_repickle or not os.path.isfile(pickle_fn):
		all_samples = get_sample_names('all')
		qp_sample_sets = calculate_qp_samples(all_samples, species_name, prev_cohort=prev_cohort)
		pickle.dump(qp_sample_sets, open(pickle_fn, 'wb'))
		return qp_sample_sets
	else:
		qp_sample_sets = pickle.load(open(pickle_fn, 'rb'))
		for cat in qp_sample_sets: # qp, non-qp, low-coverage
			old_sample_list = list(qp_sample_sets[cat])
			for sample in old_sample_list:
				if sample not in desired_samples:
					qp_sample_sets[cat].remove(sample)
		return qp_sample_sets

# ===========================================================================
# FUNCTIONS FOR DEALING WITH SAMPLE PAIRS
# ===========================================================================

# ===========================================================================
# calculate_ordered_same_subject_pairs: computes same-subject sample pairs;
# specify if timepoints should be nonconsecutive, consecutive, or longest
#
# Considers mothers and infants different subjects
#
# Returns same_subject_idxs, tuple with idx1 (lower order) and idx2 (higher)
# ===========================================================================

def calculate_ordered_same_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive'):
	idx_lower, idx_upper = [], []

	# reconstruct "timeseries" for each subject
	subject_order_idx_map = defaultdict(dict)
	for i in xrange(0,len(sample_list)):
		subject, order = sample_order_map[sample_list[i]]		 
		subject_order_idx_map[subject][order] = i
	
	# create index pairs within subjects
	for subject in subject_order_idx_map:
		sorted_orders = list(sorted(subject_order_idx_map[subject].keys()))
		
		if len(sorted_orders) < 2:	# no pairs can be formed
			continue
		
		if within_host_type=='longest':
			idx_lower.append(subject_order_idx_map[subject][sorted_orders[0]])
			idx_upper.append(subject_order_idx_map[subject][sorted_orders[-1]])
		elif within_host_type=='consecutive':
			for order_idx in xrange(1,len(sorted_orders)):
				idx_lower.append(subject_order_idx_map[subject][sorted_orders[order_idx-1]])
				idx_upper.append(subject_order_idx_map[subject][sorted_orders[order_idx]])
		elif within_host_type=='nonconsecutive':
			for order_idx_i in xrange(0,len(sorted_orders)):
				for order_idx_j in xrange(order_idx_i+1,len(sorted_orders)):
					idx_lower.append(subject_order_idx_map[subject][sorted_orders[order_idx_i]])
					idx_upper.append(subject_order_idx_map[subject][sorted_orders[order_idx_j]])
	
	same_subject_idxs = (np.array(idx_lower,dtype=np.int32), np.array(idx_upper,dtype=np.int32))		
	return same_subject_idxs

# ===========================================================================
# calculate_ordered_diff_subject_pairs: computes diff-subject sample pairs;
# only one sample considered for each subject; specify whether timepoint
# should be first or last # TODO?
#
# Returns diff_subject_idxs, tuple with idx1 and idx2
# ===========================================================================

def calculate_ordered_diff_subject_pairs(sample_order_map, sample_list=[], diff_host_type='first'):
	
	# reconstruct "timeseries" for each subject
	subject_order_idx_map = defaultdict(dict)
	for i in xrange(0,len(sample_list)):
		subject, order = sample_order_map[sample_list[i]]		 
		subject_order_idx_map[subject][order] = i
	
	sorted_subjects = sorted(subject_order_idx_map.keys()) # all subjects
	op = max if diff_host_type == 'last' else min
	idx_lower, idx_upper = [], []
	
	for subject_i_idx in xrange(0,len(sorted_subjects)):
		subject_i = sorted_subjects[subject_i_idx]
		earliest_order_i = op(subject_order_idx_map[subject_i].keys())
		i = subject_order_idx_map[subject_i][earliest_order_i]
		
		for subject_j_idx in xrange(subject_i_idx+1,len(sorted_subjects)):
			subject_j = sorted_subjects[subject_j_idx]
			earliest_order_j = op(subject_order_idx_map[subject_j].keys())
			j = subject_order_idx_map[subject_j][earliest_order_j]
			
			idx_lower.append(i)
			idx_upper.append(j)
	
	diff_subject_idxs = (np.array(idx_lower,dtype=np.int32), np.array(idx_upper,dtype=np.int32))
	return diff_subject_idxs

# ===========================================================================
# calculate_ordered_same_sample_pairs: computes same sample "pair" indices
# Assumes no duplicate samples in sample_list
# ===========================================================================

def calculate_ordered_same_sample_pairs(sample_order_map, sample_list=[]):
	idxs = np.arange(0,len(sample_list))
	return (np.array(idxs,dtype=np.int32), np.array(idxs,dtype=np.int32))

# ===========================================================================
# calculate_ordered_subject_pairs: wrapper function for combining 
# calculate_ordered_same_sample_pairs, calculate_ordered_same_subject_pairs
# and calculate_ordered_diff_subject_pairs
# ===========================================================================

def calculate_ordered_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive', diff_host_type='first'):
	
	same_sample_idxs = calculate_ordered_same_sample_pairs(sample_order_map, sample_list)
	same_subject_idxs = calculate_ordered_same_subject_pairs(sample_order_map, sample_list, within_host_type)
	diff_subject_idxs = calculate_ordered_diff_subject_pairs(sample_order_map, sample_list, diff_host_type)
	return same_sample_idxs, same_subject_idxs, diff_subject_idxs

# ===========================================================================
# calculate_mi_ordered_same_subject_pairs: computes same-subject sample pairs
# considering each mother-infant pair as the SAME subject. Parameters:
# 	within_host_type: nonconsecutive, consecutive, longest
# 	one_per_mi_pair: True, False
# 	infant_timepoint_pref: first, random, last
#
# Returns same_subject_idxs, tuple with idx1 (lower order / mother) and
# idx2 (higher order / infant)
# ===========================================================================

def calculate_mi_ordered_same_subject_pairs(sample_order_map, sample_list=[], within_host_type='consecutive', one_per_mi_pair=False, infant_timepoint_pref='first'):
		
		same_subject_idxs = calculate_ordered_same_subject_pairs(sample_order_map, sample_list, within_host_type)
		
		mother_samples = get_sample_names('mother','all')
		infant_samples = get_sample_names('infant','all')
		sample_subject_map = parse_sample_subject_map()
		
		# Include same mother and infant comparisons (one per subject)
		more_same_subject_idxs_i = []
		more_same_subject_idxs_j = []
		
		
		if one_per_mi_pair == True: # Choose one timepoint pair for each mother-infant pair
			
			existing_subjects = {}
			
			for i in range(len(sample_list)):
				for j in range(i+1, len(sample_list)):
					sample_i, sample_j = sample_list[i], sample_list[j]
					# Check if one is mother and other is her infant
					if is_mi_pair(sample_i, sample_j, mother_samples, infant_samples) and is_same_mi_subject(sample_i, sample_j, sample_subject_map):					
						m_index, i_index = (i, j) if sample_order_map[sample_i][0][-1] == 'M' else (j, i)
						subject = sample_order_map[sample_i][0][:-2] # Excludes -I/-M
						# If random, any mother-infant timepoint pair works
						if subject not in existing_subjects:
							existing_subjects[subject] = (m_index, i_index)
						elif infant_timepoint_pref != 'random':
							prev_i_index = existing_subjects[subject][1]
							prev_infant_tp = sample_order_map[sample_list[prev_i_index]][1]
							cur_infant_tp = sample_order_map[sample_list[i_index]][1]
							if infant_timepoint_pref == 'first':
								if cur_infant_tp < prev_infant_tp:
									existing_subjects[subject] = (m_index, i_index)
							elif infant_timepoint_pref == 'last':
								if cur_infant_tp > prev_infant_tp:
									existing_subjects[subject] = (m_index, i_index)
			
			for index_pair in existing_subjects.values():
				more_same_subject_idxs_i.append(index_pair[0])
				more_same_subject_idxs_j.append(index_pair[1])
		
		else: # Choose all timepoint pairs for each mother-infant pair
			
			for i in range(len(sample_list)):
				for j in range(i+1, len(sample_list)):
					sample_i, sample_j = sample_list[i], sample_list[j]
					if is_mi_pair(sample_i, sample_j, mother_samples, infant_samples) and is_same_mi_subject(sample_i, sample_j, sample_subject_map):
						m_index, i_index = (i, j) if sample_order_map[sample_i][0][-1] == 'M' else (j, i)
						more_same_subject_idxs_i.append(m_index)
						more_same_subject_idxs_j.append(i_index)
		
		idxs1 = np.array(np.append(same_subject_idxs[0], more_same_subject_idxs_i), dtype=np.int32)
		idxs2 = np.array(np.append(same_subject_idxs[1], more_same_subject_idxs_j), dtype=np.int32)
		return (idxs1, idxs2)

# ===========================================================================
# calculate_sample_idx_map: creates a map of indexes from one list of samples
# (sample_list_from) to another list of samples (sample_list_to).
# The from list must be a strict subset of the to list. 
# ===========================================================================

def calculate_sample_idx_map(sample_list_from, sample_list_to):
	
	sample_list_to = list(sample_list_to)
	sample_map = {}
	for i in xrange(0,len(sample_list_from)):
		sample_map[i] = sample_list_to.index(sample_list_from[i])
	
	return sample_map

def apply_sample_index_map_to_indices(sample_idx_map, idxs):
	new_idxs = (np.array([sample_idx_map[i] for i in idxs[0]]), np.array([sample_idx_map[i] for i in idxs[1]]))
	return new_idxs

# ===========================================================================
# calculate_samples_in_different_subjects: returns boolean array indicating
# whether each sample in sample_list has same subject as focal_sample
# ===========================================================================

def calculate_samples_in_different_subjects(sample_subject_map, sample_list, focal_sample):
	focal_subject = sample_subject_map[focal_sample]
	subjects = np.array([sample_subject_map[s] for s in sample_list])
	return (subjects != focal_subject)		

# ===========================================================================
# calculate_subject_pairs: calculates which samples belong to different 
# subjects, which belong to different timepoints in same subject, and 
# which are the same timepoint.
# 
# Considers mothers and infants different subjects
# This is a naive way to generate sample pairs
#
# Returns same_sample_idxs, same_subject_idxs, diff_subject_idxs, 
# each of which is a tuple with lists idx1 and idx2.
# All pairs are included only once (order doesn't matter).
# ===========================================================================

def calculate_subject_pairs(sample_subject_map, sample_list = None):
	
	if sample_list is None:
		sample_list = sample_subject_map.keys()
	
	same_sample_idx_lower, same_sample_idx_upper = [], []
	same_subject_idx_lower, same_subject_idx_upper = [], []
	diff_subject_idx_lower, diff_subject_idx_upper = [], []
	
	for i in xrange(0,len(sample_list)):
		sample = sample_list[i]
		same_sample_idx_lower.append(i)
		same_sample_idx_upper.append(i)
		for j in xrange(0,i):
			if sample_subject_map[sample] == sample_subject_map[sample]:
				same_subject_idx_lower.append(i)
				same_subject_idx_upper.append(j)
			else: 
				diff_subject_idx_lower.append(i)
				diff_subject_idx_upper.append(j)
		
	same_sample_idxs = (np.array(same_sample_idx_lower,dtype=np.int32), np.array(same_sample_idx_upper,dtype=np.int32))	
	same_subject_idxs = (np.array(same_subject_idx_lower,dtype=np.int32), np.array(same_subject_idx_upper,dtype=np.int32))	
	diff_subject_idxs = (np.array(diff_subject_idx_lower,dtype=np.int32), np.array(diff_subject_idx_upper,dtype=np.int32))
	
	return same_sample_idxs, same_subject_idxs, diff_subject_idxs

# ===========================================================================
# calculate_sample_subject_matrix: matrix, rows are subjects, columns are hosts 
# A_ih = 1 if sample i is in host h 
# ===========================================================================

def calculate_sample_subject_matrix(samples):
	
	sample_idx_map = {samples[i]:i for i in xrange(0,len(samples))}
	
	subject_sample_map = parse_subject_sample_map()
	subjects = subject_sample_map.keys()
	
	sample_subject_matrix = np.zeros((len(samples),len(subjects)),dtype=np.bool)
	
	for subject_idx in xrange(0,len(subjects)):
		for sample in subject_sample_map[subjects[subject_idx]]:
			if sample in sample_idx_map:
				sample_subject_matrix[sample_idx_map[sample], subject_idx] = True
	
	return sample_subject_matrix, subjects

# ===========================================================================
# FUNCTIONS FOR SAMPLE-METADATA MAPS
# ===========================================================================

# ====================================================================================
# parse_sample_metadata_map
# 
# Loads metadata for HMP, Yassour, Backhed, Ferretti, Shao, Olm samples
# Returns map:
# sample -> (subject_id, sample_id, accession_id, country, continent, temporal_order)
# 
# By default (at least for postprocessing steps), include all samples. But can
# restrict to fecal samples only, and good timepoints only (some Shao have 'NA' tp)
# ====================================================================================

def parse_sample_metadata_map(fecal_only = False, good_tp_only = False): 
	
	metadata_dir = config.metadata_directory
	samples_dir = "%s/final_sample_lists" % metadata_dir
	
	sample_metadata_map = {} # What is returned!
	
	# First load HMP metadata (469 samples)
	
	samples_fpath = "%s/HMP1-2_samples.txt" % samples_dir
	hmp_samples = [line.strip() for line in open(samples_fpath, 'r')]
	hmp_samples = parse_merged_sample_names(hmp_samples) # Remove c's
	
	with open("%s/HMP1-2_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			subject_id, sample_id, accession_id, country, continent, order = line.strip().split('\t')
			order = int(order)
			if sample_id in hmp_samples:
				sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)
	
	# Then load Backhed data (391 samples)
	
	timept_order_map_mother = {'M':1}
	timept_order_map_infant = {'B':1,'4M':2,'12M':3}
	
	with open("%s/Backhed_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			accession_id, subject_id, timept = line.strip().split('\t')
			# Using the family/study_id as subject id, and specify mother (-M) vs infant (-I)
			subject_id = subject_id + ('-M' if timept == 'M' else '-I')
			sample_id = accession_id # Sample ID same as run accession
			order = timept_order_map_mother[timept] if timept == 'M' else timept_order_map_infant[timept]
			sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Sweden', 'Europe', order)
	
	# Then load Ferretti data (119 fecal samples)
	
	timept_order_map_mother = {'t0':1}
	timept_order_map_infant = {'t1':1,'t2':2,'t3':3,'t4':4,'t5':5}
	
	with open("%s/Ferretti_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			items = line.strip().split("\t")
			accession_id = items[4]
			subject_id = items[6][:11] # Using first 11 characters of experiment_alias as subject id (e.g. CA_C10055IS)
			sample_id = accession_id # Sample ID same as run accession
			timept = items[5][-2:] # Using last two characters of experiment_alias to identify timepoint
			order = timept_order_map_mother[timept] if subject_id[-2:] == 'MS' else timept_order_map_infant[timept]
			if fecal_only:
				if items[6][15:17] == 'FE': # Restrict to fecal samples
					sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Italy', 'Europe', order)
			else:
				sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Italy', 'Europe', order)
	
	# Then load Yassour data (286 samples)
	
	timept_order_map_mother = {'Mother:Gest':1, 'Mother:Birth':2,'Mother:3 months':3}
	timept_order_map_infant = {'Child:Birth':1,'Child:14 days':2,'Child:1 month':3,'Child:2 months':4,'Child:3 months':5}
	
	with open("%s/Yassour_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			items = line.strip().split("\t")
			accession_id = items[4]
			subject_id = items[7][:7] # Using first 7 characters of sample_title as subject id (e.g. M0059-M)
			sample_id = accession_id # Sample ID same as run accession
			timept = items[7][6:] # Using characters after 6th of sample_title to identify timepoint
			order = timept_order_map_mother[timept] if subject_id[-1] == 'M' else timept_order_map_infant[timept]
			sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'Finland', 'Europe', order)
	
	# Then load Shao data (1679 samples - 1676 excluding NA timepoints)
	
	with open("%s/Shao_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			accession_id, _, subject_id, timept, infancy_months = line.strip().split("\t")
			sample_id = accession_id # Sample ID same as run accession
			
			# Adjust subject ID by appending -M or -I
			if timept == 'Mother' and infancy_months == 'Mother':
				subject_id += '-M'
				order = 1 # Only one mother timepoint
			elif timept == 'Infancy':
				subject_id += '-I'
				
				if infancy_months == 'NA':
					if good_tp_only == True:
						continue # Skip if infancy months field is NA
					else:
						infancy_months = -999 # Bogus negative number
				
				order = float(infancy_months) * 30.5 # Convert months to approx. days
			elif infancy_months == 'Neonatal':
				subject_id += '-I'
				order = int(timept) # In days since birth
			
			sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'United Kingdom', 'Europe', order)
	
	# Then load Olm data (898 samples)
	
	olm_samples = []
	for campaign in ['NIH1', 'NIH2', 'NIH3', 'NIH4', 'Sloan2']:
		samples_fpath = "%s/Olm_%s_samples.txt" % (samples_dir, campaign)
		olm_sub_samples = [line.strip() for line in open(samples_fpath, 'r')]
		olm_samples += list(parse_merged_sample_names(olm_sub_samples)) # Remove c's
	
	with open("%s/Olm_metadata.txt" % metadata_dir, 'r') as metadata_file:
		metadata_file.readline() # header
		for line in metadata_file:
			items = line.strip().split("\t")
			if len(items) == 10: # Must have available accession
				subject_id = items[1]
				timept = items[2]
				accession_id = items[9]
				sample_id = accession_id # Sample ID same as run accession
				order = int(timept) # In days since birth
				if accession_id in olm_samples: # Restrict to considered samples
					sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, 'United States', 'North America', order)
	
	return sample_metadata_map

# ====================================================================================
# filter_sample_metadata_map
# 
# Using passed in sample-metadata map, filters only sample entries corresponding to a
# certain subject_id, country, continent or order
# ====================================================================================

def filter_sample_metadata_map(sample_metadata_map, field, field_value):
	
	field_dict = {"subject_id": 0, "country": 3, "continent": 4, "order": 5}
	if field in field_dict:
		field_idx = field_dict[field]
	else:
		return sample_metadata_map
	
	filtered_sample_metadata_map = {}
	for sample in sample_metadata_map:
		if sample_metadata_map[sample][field_idx] == field_value:
			filtered_sample_metadata_map[sample] = sample_metadata_map[sample]
	
	return filtered_sample_metadata_map

# ====================================================================================
# extract_sample_metadata_map
# 
# Using passed in sample-metadata map, extracts information for only one column
# (options: subject_id, country, continent or order) and returns new map
# Loads the default full sample-metadata map if nothing passed in
# ====================================================================================

def extract_sample_metadata_map(field, sample_metadata_map = None):
	
	field_dict = {"subject_id": 0, "country": 3, "continent": 4, "order": 5}
	field_idx = field_dict[field] if (field in field_dict) else 0 # Defaults to subject_id
	
	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it
	
	extracted_sample_metadata_map = {}
	for sample in sample_metadata_map:
		extracted_sample_metadata_map[sample] = sample_metadata_map[sample][field_idx]
	
	return extracted_sample_metadata_map

# ====================================================================================
# parse_sample_subject_map, parse_sample_country_map, parse_sample_continent_map
#
# Convenience functions for extract_sample_metadata_map
# ====================================================================================

def parse_sample_subject_map(sample_metadata_map = None): 
	return extract_sample_metadata_map("subject_id", sample_metadata_map)

def parse_sample_country_map(sample_metadata_map = None): 
	return extract_sample_metadata_map("country", sample_metadata_map)

def parse_sample_continent_map(sample_metadata_map = None): 
	return extract_sample_metadata_map("continent", sample_metadata_map)

# ====================================================================================
# parse_subject_sample_map
# 
# Returns map: subject -> map: samples -> set of accession IDs
# ====================================================================================

def parse_subject_sample_map(sample_metadata_map = None): 
		
	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it
	
	subject_sample_map = {}
	for sample in sample_metadata_map:
		subject_id, _, accession_id, country, continent, order = sample_metadata_map[sample]		
		if subject_id not in subject_sample_map:
			subject_sample_map[subject_id] = {}				
		if sample not in subject_sample_map[subject_id]:
			subject_sample_map[subject_id][sample] = set()		
		subject_sample_map[subject_id][sample].add(accession_id)
	
	return subject_sample_map

# ====================================================================================
# parse_sample_order_map
# 
# Returns map from sample -> (subject_id, temporal_order)
# ====================================================================================

def parse_sample_order_map(sample_metadata_map = None): 
	
	if sample_metadata_map is None:
		sample_metadata_map = parse_sample_metadata_map() # Load it
	
	sample_order_map = {}
	for sample in sample_metadata_map:
			subject_id, _, _, _, _, order = sample_metadata_map[sample]
			sample_order_map[sample] = (subject_id, order)
	
	return sample_order_map

# ====================================================================================
# parse_sample_cohort_map
# ====================================================================================

def parse_sample_cohort_map():
	
	sample_cohort_map = {}
	sample_dict = get_sample_names('all-dict')
	
	for cohort in sample_dict:
		for sample_list in sample_dict[cohort].values():
			for sample in sample_list:
				sample_cohort_map[sample] = cohort
	
	return sample_cohort_map

# ====================================================================================
# FUNCTIONS FOR DEALING WITH REPLICATES		
# ====================================================================================

# ====================================================================================
# Returns a flat map of all the replicate sets for
# the samples in subject_sample_map, indexed by sample key				
# ====================================================================================

def flatten_samples(subject_sample_map):
	
	grouping_replicate_map = {}
	for subject in sorted(subject_sample_map.keys()):
		for sample in sorted(subject_sample_map[subject].keys()):
			grouping_replicate_map[sample] = subject_sample_map[subject][sample]
	
	return grouping_replicate_map

# ====================================================================================
# Returns a flat map of the merged replicate sets for each subject, 
# indexed by subject key 
# ====================================================================================
 
def flatten_subjects(subject_sample_map):
	
	grouping_replicate_map = {}
	for subject in sorted(subject_sample_map.keys()):
		merged_replicates = set()
		for sample in subject_sample_map[subject].keys():
			merged_replicates.update(subject_sample_map[subject][sample])
		grouping_replicate_map[subject] = merged_replicates
	
	return grouping_replicate_map

# ====================================================================================
# groupings = ordered list of nonoverlapping sets of sample names
# samples = ordered list of samples
#
# returns: list whose i-th element contains a np array of idxs
#					 of the items in samples that are present in the ith grouping
# ====================================================================================
			
def calculate_grouping_idxs(groupings, samples):
		
		grouping_idxs = []
		for i in xrange(0,len(groupings)):
		
				idxs = []
				for j in xrange(0,len(samples)):
						if samples[j] in groupings[i]:
								idxs.append(j)
				idxs = np.array(idxs,dtype=np.int32)
				#print idxs
				grouping_idxs.append(idxs)
		
		return grouping_idxs

# ===========================================================================
# Isolate and mixture metadata parsing
# ===========================================================================

def parse_isolate_metadata_map():
		
		isolate_metadata_map = {}
		
		# load simulations
		file = open(parse_midas_data.scripts_directory+"isolates_genome_list.txt","r")
		file.readline() # 
		for line in file:
				items = line.strip().split("\t")
				subject_id = items[0] 
				sample_id = subject_id 
				accession_id=subject_id
				country = "isolate"
				continent = "isolate"
				order = 1
				
				isolate_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)
				
		file = open(parse_midas_data.scripts_directory+"mixture_labels.txt","r")
		file.readline() # header
		for line in file:
				items = line.strip().split("\t")
				subject_id = items[0] # this is one of two of the 90/10 mixtures
				sample_id = items[1] # This is the exact simulation
				accession_id=sample_id # same as sample
				country = "mixture"
				continent = "mixture"
				order = 1
				isolate_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)
				
		return isolate_metadata_map

def list_of_isolates_and_mixtures():
		
		isolate_metadata_map = parse_isolate_metadata_map()
		
		isolates=[]
		mixtures=[]
		
		for sample_id in isolate_metadata_map:
				subject_id, dummy, accession_id, country, continent, order = sample_metadata_map[sample_id]
				if country=='isolate':
						isolates.append(sample_id)
				elif country=='mixture':
						mixtures.append(sample_id)
						 
		return isolates, mixtures

# ===========================================================================
# Simply removes 'c' suffix for any merged samples
# ===========================================================================

def parse_merged_sample_names(items):
	samples = []
	for item in items:
		sample = item.strip()
		if sample.endswith('c'):
			sample = sample[:-1]
		samples.append(sample)
	
	samples = np.array(samples)
	return samples

# ===========================================================================
# Prunes sample list to remove multiple timepoints from same subject
# Considers mothers and infants as different subjects
# Returns len(sample_list) boolean array with element=False if sample was pruned	
# ===========================================================================

def calculate_unique_samples(subject_sample_map, sample_list=None):

		if sample_list is None:
				sample_list = list(sorted(flatten_samples(subject_sample_map).keys()))
		
		# invert subject sample map
		sample_subject_map = parse_sample_subject_map()
		
		subject_idx_map = {}
		
		for i in xrange(0,len(sample_list)):
				sample = sample_list[i]
				if sample.endswith('c'):
						sample = sample[:-1]
				subject = sample_subject_map[sample]
				if not subject in subject_idx_map:
						subject_idx_map[subject] = i
						
		unique_idxs = np.zeros(len(sample_list),dtype=np.bool_)
		for i in subject_idx_map.values():
				unique_idxs[i]=True
		
		return unique_idxs
