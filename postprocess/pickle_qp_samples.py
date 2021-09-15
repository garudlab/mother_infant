import pickle
import sample_utils as su
import parse_midas_data

good_species_list = parse_midas_data.parse_good_species_list()
cohorts = ['hmp', 'backhed', 'yassour', 'ferretti']

sample_country_map = su.parse_sample_country_map()
all_samples = sample_country_map.keys()
hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

# Stores indices with respect to all_samples
qp_samples = {cohort: {} for cohort in cohorts}

for species_name in good_species_list:
	print(species_name)
	for cohort in cohorts:
		if cohort == 'hmp':
			this_qp_samples = su.calculate_qp_samples(hmp_samples, species_name)['qp']
		else:
			this_qp_samples = su.calculate_qp_samples(su.get_sample_names(cohort,'all'), species_name)['qp']
		qp_indices = [all_samples.index(s) for s in this_qp_samples]
		qp_samples[cohort][species_name] = qp_indices

import config

all_samples_pickle_fn = '%s/pickles/all_samples.pkl' % config.data_directory
pickle.dump(all_samples, open(all_samples_pickle_fn, 'w'))

qp_samples_pickle_fn = '%s/pickles/qp_samples.pkl' % config.data_directory
pickle.dump(qp_samples, open(qp_samples_pickle_fn, 'w'))
