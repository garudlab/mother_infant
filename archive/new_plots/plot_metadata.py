from utils import sample_utils as su, config
import numpy as np
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

sample_order_map = su.parse_sample_order_map()
hmp_samples = su.get_sample_names('hmp')
mother_samples = su.get_sample_names('mother')

# ==============================================================
# Backhed timepoints:
# M: Mother, delivery / 0-5 days after birth (median 2)
# B: Infant, birth / 2-5 days after birth (median 3)
# 12M: Infant, 4 months / 119-125 days after birth (median 122)
# 4M: Infant, 12 months / 363-372 days after birth (median 366)
# ==============================================================

cohort = 'Backhed'
subject_tps_map = defaultdict(list) # subject -> list of timepoints

for sample in su.get_sample_names('backhed'):
	subject, order = sample_order_map[sample]
	tp = su.sample_to_tp(sample, sample_order_map, hmp_samples, mother_samples)
	subject_tps_map[subject].append(tp)

subject_tps_map = {subject: sorted(orders) for subject, orders in subject_tps_map.items()}
subject_tps = sorted(subject_tps_map.items(), key=lambda x: x[1])
subjects, all_tps = zip(*subject_tps)

subject_idx = 0
subject_idxs = []
flattened_tps = []
for tps in all_tps:
	for tp in tps:
		subject_idxs.append(subject_idx)
		flattened_tps.append(tp)
	subject_idx += 1

tp_idx_map = {'I1': 1, 'I2': 2, 'I3': 3, 'M1': 0}
tp_label_map = {'I1': '~3 days', 'I2': '4 months', 'I3': '12 months', 'M1': 'Mother\n~2 days'}

flattened_tp_idxs = [tp_idx_map[tp] for tp in flattened_tps]

# Try points...?
fig, ax = plt.subplots(figsize=(6, 7))
ax.plot(flattened_tp_idxs, subject_idxs, '_', markersize=50)
ax.set_xlim((-0.5, 3.5))
ax.set_xticks((0, 1, 2, 3))
ax.set_xticklabels([tp_label_map[tp] for tp in ('M1', 'I1', 'I2', 'I3')])
ax.set_title('Backhed timepoints\n(%i samples, %i subjects)' % (len(subject_idxs), len(set(subject_idxs))))

fig.savefig("%s/%s_metadata.pdf" % (config.analysis_directory, cohort), bbox_inches='tight')

'''
# Try heatmap...?
heatmap = np.zeros((len(tp_idx_map), len(subject_idxs)))
for t_idx, s_idx in zip(flattened_tp_idxs, subject_idxs):
	heatmap[t_idx][s_idx] = (t_idx + 1)
	
plt.imshow(heatmap)
plt.show()
'''
# ==============================================================
# Ferretti timepoints:
# M0: Mother, delivery
# I1: Infant, 1 day
# I2: Infant, 3 days
# I3: Infant, 1 week
# I4: Infant, 1 month
# I5: Infant, 4 months
# ==============================================================

cohort = 'Ferretti'
subject_tps_map = defaultdict(list) # subject -> list of timepoints

for sample in su.get_sample_names(cohort):
	subject, order = sample_order_map[sample]
	tp = su.sample_to_tp(sample, sample_order_map, hmp_samples, mother_samples)
	subject_tps_map[subject].append(tp)

subject_tps_map = {subject: sorted(orders) for subject, orders in subject_tps_map.items()}
subject_tps = sorted(subject_tps_map.items(), key=lambda x: x[1])
subjects, all_tps = zip(*subject_tps)

subject_idx = 0
subject_idxs = []
flattened_tps = []
for tps in all_tps:
	for tp in tps:
		subject_idxs.append(subject_idx)
		flattened_tps.append(tp)
	subject_idx += 1

tp_idx_map = {'I1': 1, 'I2': 2, 'I3': 3, 'I4': 4, 'I5': 5, 'M1': 0}
tp_label_map = {'I1': '1 day', 'I2': '3 days', 'I3': '1 week', 'I4': '1 month', 'I5': '4 months', 'M1': 'Mother\ndelivery'}

flattened_tp_idxs = [tp_idx_map[tp] for tp in flattened_tps]

# Plot
fig, ax = plt.subplots(figsize=(6, 7))
ax.plot(flattened_tp_idxs, subject_idxs, '_', markersize=50)
ax.set_xlim((-0.5, 5.5))
ax.set_xticks((0, 1, 2, 3, 4, 5))
ax.set_xticklabels([tp_label_map[tp] for tp in ('M1', 'I1', 'I2', 'I3', 'I4', 'I5')])
ax.set_title('Ferretti timepoints\n(%i samples, %i subjects)' % (len(subject_idxs), len(set(subject_idxs))))

fig.savefig("%s/%s_metadata.pdf" % (config.analysis_directory, cohort), bbox_inches='tight')

# ==============================================================
# Yassour timepoints:
# MGest: Mother, gestational week 27
# MBirth: Mother, delivery
# M3: Mother, 3 months post-delivery
# CBirth: Infant, birth / meconium
# C14: Infant, 2 weeks
# C1: Infant, 1 month
# C2: Infant, 2 months
# C3: Infant, 3 months
# ==============================================================

cohort = 'Yassour'
subject_tps_map = defaultdict(list) # subject -> list of timepoints

for sample in su.get_sample_names(cohort):
	subject, order = sample_order_map[sample]
	tp = su.sample_to_tp(sample, sample_order_map, hmp_samples, mother_samples)
	subject_tps_map[subject].append(tp)

subject_tps_map = {subject: sorted(orders) for subject, orders in subject_tps_map.items()}
subject_tps = sorted(subject_tps_map.items(), key=lambda x: x[1])
subjects, all_tps = zip(*subject_tps)

subject_idx = 0
subject_idxs = []
flattened_tps = []
for tps in all_tps:
	for tp in tps:
		subject_idxs.append(subject_idx)
		flattened_tps.append(tp)
	subject_idx += 1

tp_idx_map = {'I1': 3, 'I2': 4, 'I3': 5, 'I4': 6, 'I5': 7, 'M1': 0, 'M2': 1, 'M3': 2}
tp_label_map = {'I1': 'Meconium', 'I2': '2 weeks', 'I3': '1 month', 'I4': '2 months', 'I5': '3 months', 'M1': 'Mother\nGest wk27', 'M2': 'Mother\ndelivery', 'M3': 'Mother\n3 months'}

flattened_tp_idxs = [tp_idx_map[tp] for tp in flattened_tps]

# Plot
fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(flattened_tp_idxs, subject_idxs, '_', markersize=50)
ax.set_xlim((-0.5, 7.5))
ax.set_xticks((0, 1, 2, 3, 4, 5, 6, 7))
ax.set_xticklabels([tp_label_map[tp] for tp in ('M1', 'M2', 'M3', 'I1', 'I2', 'I3', 'I4', 'I5')])
ax.set_title('Yassour timepoints\n(%i samples, %i subjects)' % (len(subject_idxs), len(set(subject_idxs))))

fig.savefig("%s/%s_metadata.pdf" % (config.analysis_directory, cohort), bbox_inches='tight')

# ==============================================================
# Shao timepoints:
# Mother: one timepoint at delivery
# Neonatal: 4-21 days (various)
# Infancy: 4-15 months (various)
# ==============================================================

cohort = 'Shao'
subject_tps_map = defaultdict(list) # subject -> list of timepoints

for sample in su.get_sample_names(cohort):
	subject, order = sample_order_map[sample]
	tp = su.sample_to_tp(sample, sample_order_map, hmp_samples, mother_samples)
	subject_tps_map[subject].append(tp)

subject_tps_map = {subject: sorted(orders) for subject, orders in subject_tps_map.items()}
subject_tps = sorted(subject_tps_map.items(), key=lambda x: x[1])
subjects, all_tps = zip(*subject_tps)

subject_idx = 0
subject_idxs = []
flattened_tps = []
for tps in all_tps:
	for tp in tps:
		subject_idxs.append(subject_idx)
		flattened_tps.append(tp)
	subject_idx += 1

tp_idx_map = {tp: int(np.round(float(tp[1:]))) for tp in set(flattened_tps)}
sorted_tp_idxs = sorted(tp_idx_map.values())

flattened_tp_idxs = [tp_idx_map[tp] for tp in flattened_tps]

# Plot
fig, ax = plt.subplots(figsize=(15, 24))
ax.plot(flattened_tp_idxs, subject_idxs, '_', markersize=5)
ax.set_xlim((-0.5, max(sorted_tp_idxs) + 0.5))
ax.set_xlabel("Days (first tick is mother, delivery)")
# ax.set_xticks(sorted_tp_idxs)
# ax.set_xticklabels(['M'] + [str(tp) for tp in sorted_tp_idxs])
ax.set_title('Shao timepoints\n(%i samples, %i subjects)' % (len(subject_idxs), len(set(subject_idxs))))

fig.savefig("%s/%s_metadata.pdf" % (config.analysis_directory, cohort), bbox_inches='tight')

# ==============================================================
# Olm timepoints: 5-86 days
# ==============================================================

cohort = 'Olm'
subject_tps_map = defaultdict(list) # subject -> list of timepoints

for sample in su.get_sample_names(cohort):
	subject, order = sample_order_map[sample]
	tp = su.sample_to_tp(sample, sample_order_map, hmp_samples, mother_samples)
	subject_tps_map[subject].append(tp)

subject_tps_map = {subject: sorted(orders) for subject, orders in subject_tps_map.items()}
subject_tps = sorted(subject_tps_map.items(), key=lambda x: x[1])
subjects, all_tps = zip(*subject_tps)

subject_idx = 0
subject_idxs = []
flattened_tps = []
for tps in all_tps:
	for tp in tps:
		subject_idxs.append(subject_idx)
		flattened_tps.append(tp)
	subject_idx += 1

tp_idx_map = {tp: int(tp[1:]) for tp in set(flattened_tps)}
sorted_tp_idxs = sorted(tp_idx_map.values())

flattened_tp_idxs = [tp_idx_map[tp] for tp in flattened_tps]

# Plot
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(flattened_tp_idxs, subject_idxs, '_', markersize=5)
ax.set_xlim((-0.5, max(sorted_tp_idxs) + 0.5))
ax.set_xlabel("Days")
ax.set_title('Olm timepoints\n(%i samples, %i subjects)' % (len(subject_idxs), len(set(subject_idxs))))

fig.savefig("%s/%s_metadata.pdf" % (config.analysis_directory, cohort), bbox_inches='tight')
