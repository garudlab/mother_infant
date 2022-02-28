from utils import parse_midas_data, sample_utils, config
from collections import defaultdict

# Load species list, timepoint-sample dict
good_species_list = parse_midas_data.load_pickled_good_species_list()
mi_tp_sample_dict = sample_utils.get_mi_tp_sample_dict(exclude_cohorts = ['olm'])
infant_tps_ordered = sorted(mi_tp_sample_dict['infant'].keys())

# =======================================================================
# Build qp information for infants (species -> tp -> qp samples count)
# =======================================================================

num_qp_dict = defaultdict(dict)
num_non_dict = defaultdict(dict)
num_lowcov_dict = defaultdict(dict)

for species in good_species_list:
	
	print("Working on species %s..." % species)
	
	for tp in mi_tp_sample_dict['infant']:
		
		samples = mi_tp_sample_dict['infant'][tp]
		num_samples = len(samples)
		
		qp_sample_sets = sample_utils.load_qp_samples(samples, species)
		num_qp = len(qp_sample_sets['qp'])
		num_non = len(qp_sample_sets['non-qp'])
		num_lowcov = len(qp_sample_sets['low-coverage'])
		
		num_qp_dict[species][tp] = num_qp
		num_non_dict[species][tp] = num_non
		num_lowcov_dict[species][tp] = num_lowcov

# =======================================================================
# Pickle
# =======================================================================

import pickle
ddir = config.data_directory

pickle.dump(num_qp_dict, open("%s/pickles/plot_qp/num_qp_dict.pkl" % ddir, 'wb'))
pickle.dump(num_non_dict, open("%s/pickles/plot_qp/num_non_dict.pkl" % ddir, 'wb'))
pickle.dump(num_lowcov_dict, open("%s/pickles/plot_qp/num_lowcov_dict.pkl" % ddir, 'wb'))

# =======================================================================
# Bar plots comparing QP vs. non-QP sample counts over timepoints
# aggregated over species + combining datasets
# =======================================================================

import numpy as np
from matplotlib import pyplot as plt

num_qp_agg_species = []
num_non_agg_species = []
num_lowcov_agg_species = []

for tp in infant_tps_ordered:
	total_num_qp = 0
	total_num_non = 0
	total_num_lowcov = 0
	
	for species in good_species_list:
		total_num_qp += num_qp_dict[species][tp]
		total_num_non += num_non_dict[species][tp]
		total_num_lowcov += num_lowcov_dict[species][tp]
	
	num_qp_agg_species.append(total_num_qp)
	num_non_agg_species.append(total_num_non)
	num_lowcov_agg_species.append(total_num_lowcov)

labels = [str(tp) for tp in infant_tps_ordered]
xticks = np.arange(len(labels))

fig, ax = plt.subplots(figsize=(28,8))
ax.set_yscale('log')
ax.set_xlim((-0.5, max(infant_tps_ordered) + 0.5))
ax.bar(infant_tps_ordered, num_qp_agg_species, label='QP', color='orange', linewidth=0)
ax.bar(infant_tps_ordered, num_non_agg_species, label='non-QP', bottom=num_qp_agg_species, color='blue', linewidth=0)
num_highcov_agg_species = np.array(num_qp_agg_species) + np.array(num_non_agg_species)
ax.bar(infant_tps_ordered, num_lowcov_agg_species, label='low cov', bottom=num_highcov_agg_species, color='gray', linewidth=0)
# ax.set_xticklabels(labels)
ax.hlines(10, -0.5, max(infant_tps_ordered) + 0.5), linestyles='dashed')
ax.set_ylabel("Number of samples (from any cohort)")
ax.set_xlabel("Timepoint (days)")
ax.set_title("Number of QP samples by timepoint")
ax.legend()

fig.savefig("%s/num_qp_over_time.pdf" % config.analysis_directory, bbox_inches='tight')

# =======================================================================
# Bar/line plots comparing QP vs. non-QP proportions over timepoints
# =======================================================================

prop_qp_agg_species = []
prop_non_agg_species = []
new_infant_tps_ordered = []

for q, n, tp in zip(num_qp_agg_species, num_non_agg_species, infant_tps_ordered):
	num_total = float(q + n)
	if num_total != 0:
		prop_qp_agg_species.append(q/num_total)
		prop_non_agg_species.append(n/num_total)
		new_infant_tps_ordered.append(tp)

fig, ax = plt.subplots(figsize=(28,8))
ax.set_xlim((-0.5, max(new_infant_tps_ordered) + 0.5))
ax.bar(new_infant_tps_ordered, prop_qp_agg_species, label='QP', color='orange', linewidth=0)
ax.bar(new_infant_tps_ordered, prop_non_agg_species, label='non-QP', bottom=prop_qp_agg_species, color='blue', linewidth=0)
ax.set_ylabel("Proportion of high coverage samples (from any cohort)\n")
ax.set_xlabel("Timepoint (days)")
ax.set_title("Proportion of QP samples by timepoint")
ax.legend()

fig.savefig("%s/prop_qp_over_time.pdf" % config.analysis_directory, bbox_inches='tight')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(new_infant_tps_ordered, prop_qp_agg_species)
ax.set_ylabel("Proportion of high coverage samples with are QP (from any cohort)\n")
ax.set_xlabel("Timepoint (days)")
ax.set_title("Proportion of QP samples by timepoint")
fig.savefig("%s/prop_qp_over_time_line.pdf" % config.analysis_directory, bbox_inches='tight')

# =======================================================================
# Line plots comparing QP vs. non-QP proportion over timepoints
# for a limited subset of species
# =======================================================================

prop_qp_by_species = defaultdict(list)
num_highcov_by_species = defaultdict(list)
	
for species in good_species_list[0:10]:
	
	for tp in infant_tps_ordered:
		num_qp, num_non = num_qp_dict[species][tp], num_non_dict[species][tp]
		total_highcov = num_qp + num_non
		
		if total_highcov == 0:
			prop_qp_by_species[species].append(None)
		else:
			prop_qp = num_qp / float(total_highcov)
			prop_qp_by_species[species].append(prop_qp)
		
		num_highcov_by_species[species].append(total_highcov)
