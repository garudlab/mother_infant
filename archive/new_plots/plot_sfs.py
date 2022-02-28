from collections import defaultdict
import parse_midas_data, sample_utils as su
import numpy

good_species_list =  parse_midas_data.parse_good_species_list()
species_name = good_species_list[0]
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D']))

mother_samples = su.get_sample_names('mother', 'all')
infant_samples = su.get_sample_names('infant', 'all')
hmp_samples = su.get_sample_names('hmp', 'all')

cohort_dict = {'mother': mother_samples, 'infant': infant_samples, 'adult': hmp_samples}

freq_counts = {cohort: defaultdict(int) for cohort in cohort_dict.keys()}

def get_sample_cohort(sample):
	for cohort in cohort_dict:
		if sample in cohort_dict[cohort]:
			return cohort
	
	print("Error finding cohort for sample " + str(sample))
	return None

for sample in sfs_map:
	cohort = get_sample_cohort(sample)
	if cohort is None:
		continue
	for depth, alt in sfs_map[sample]:
		f = alt*1.0/depth
		n, _ = sfs_map[sample][(depth, alt)]
		freq_counts[cohort][f] += n

binned_freq_list = np.zeros(50).astype('int')

for cohort in cohort_dict:
	for freq in sorted(freq_counts.keys()):
		if freq > 0.5:
			freq = 1 - freq
		idx = int(np.round(freq*100))
		binned_freq_list[idx] += freq_counts[freq]