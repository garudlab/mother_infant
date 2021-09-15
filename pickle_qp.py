from utils import parse_midas_data, sample_utils as su, config
from collections import defaultdict

# "Categories" (may include Olm, HMP later)
cats = ['mother', 'infant', 'hmp']

# Load species and samples list
good_species_list = parse_midas_data.load_pickled_good_species_list()
samples = {cat: su.get_sample_names(cat) for cat in cats}

# =======================================================================
# Consolidate QP information
# =======================================================================

# Map: mother / infant -> samples
sample_species_qp_dict = {cat: defaultdict(dict) for cat in cats}

for cat in cats:
	for species in good_species_list:
		print("Working on species %s..." % species)	
		qp_sample_sets = su.load_qp_samples(samples[cat], species)
		for qp_status in ['qp', 'non-qp', 'low-coverage']:
			for sample in qp_sample_sets[qp_status]:
				sample_species_qp_dict[cat][sample][species] = qp_status			

# =======================================================================
# Pickle
# =======================================================================

import pickle
pickle.dump(sample_species_qp_dict, open("%s/pickles/sample_species_qp_dict.pkl" % config.data_directory, 'wb'))
