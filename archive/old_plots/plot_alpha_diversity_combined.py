# Backhed metadata: /u/project/ngarud/daisyche/scripts/MotherInfant.run_to_sample.txt
# Run accessions all start with 'ERR'

# Yassour metadata: /u/project/ngarud/daisyche/scripts/PRJNA475246.txt
# Run accessions all start with 'SRR728'

# Ferretti metadata: /u/project/ngarud/daisyche/scripts/PRJNA352475.txt
# Run accessions all start with 'SRR527'

import sys
import bz2
import math
import config
# import numpy as np
import matplotlib	 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
		
data_fpath = '%s/species/relative_abundance.txt.bz2' % config.data_directory
cohorts = ['backhed', 'yassour', 'ferretti']
metadata_fpath = {'backhed': '%s/MotherInfant.run_to_sample.txt' % config.metadata_directory, 
'ferretti': '%s/PRJNA352475.txt' % config.metadata_directory, 
'yassour': '%s/PRJNA475246.txt' % config.metadata_directory}

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(16,4), gridspec_kw={'width_ratios': [1.4, 2, 2.4], 'wspace': 0, 'hspace': 0})

# Backhed
sys.stderr.write("Computing alpha diversity for Backhed\n")

with open(data_fpath, 'r') as data_file:
	decompressor = bz2.BZ2Decompressor()
	raw = decompressor.decompress(data_file.read())
data = [row.split('\t') for row in raw.split('\n')]
data.pop() # Get rid of extra element due to terminal newline
header = data[0]
with open(metadata_fpath['backhed'], 'r') as metadata_file:
	metadata = [row.strip().split('\t') for row in metadata_file.readlines()]
	
shannon_div_dict = {}
species_num_dict = defaultdict(int)
for i in range(len(header)):
	sample = header[i]
	if sample.startswith("ERR"):
		acc = 0
		for row in data[1:]:
			try:
				rel_ab = float(row[i])
			except ValueError:
				sys.stderr.write('Non-numeric data error\n')
				sys.exit(1)
			if rel_ab != 0:
				acc += (rel_ab * math.log(rel_ab))
				species_num_dict[sample] += 1
		shannon_div_dict[sample] = (acc*-1)

# Unordered lists of all samples' alpha diversities in given cohort
birth_alpha_div = []
month4_alpha_div = []
month12_alpha_div = []
mother_alpha_div = []

for sample in metadata[1:]:
	if sample[2] == 'M':
		mother_alpha_div.append(species_num_dict[sample[0]])
	elif sample[2] == 'B':
		birth_alpha_div.append(species_num_dict[sample[0]])
	elif sample[2] == '4M':
		month4_alpha_div.append(species_num_dict[sample[0]])
	elif sample[2] == '12M':
		month12_alpha_div.append(species_num_dict[sample[0]])
		
plot_data = [mother_alpha_div, birth_alpha_div, month4_alpha_div, month12_alpha_div]
plot_names = ["Mothers\nDelivery", "3-5 Days","4 Months", "12 Months"]

ax[0].boxplot(plot_data)
# ax[0].set_ylabel("Shannon Diversity Index")
ax[0].set_ylabel("Number of species")
ax[0].set_xticklabels(plot_names, fontsize=9)
ax[0].set_title("Backhed")
# plt.legend([bp1["boxes"][i] for i in range(0,4)], plot_names, loc='upper left')

# Ferretti
sys.stderr.write("Computing alpha diversity for Ferretti\n")

with open(metadata_fpath['ferretti'], 'r') as metadata_file:
	metadata = [row.strip().split('\t') for row in metadata_file.readlines()]
	
shannon_div_dict = {}
species_num_dict = defaultdict(int)
for i in range(len(header)):
	sample = header[i]
	if sample.startswith("SRR527"):
		acc = 0
		for row in data[1:]:
			try:
				rel_ab = float(row[i])
			except ValueError:
				sys.stderr.write('Non-numeric data error\n')
				sys.exit(1)
			if rel_ab != 0:
				acc += (rel_ab * math.log(rel_ab))
				species_num_dict[sample] += 1
		shannon_div_dict[sample] = (acc*-1)

# Unordered lists of all samples' alpha diversities in given cohort
# t0 = delivery, t1 = 1 day after, t2 = 3 days after, t3 = 1 week after, t4 = 1 month after, t3 = 4 months after
mom_t0_alpha_div = []	
child_t1_alpha_div = []
child_t2_alpha_div = []
child_t3_alpha_div = []
child_t4_alpha_div = []
child_t5_alpha_div = []

for sample in metadata[1:]:
	if sample[6][9:11] == 'MS' and sample[6][15:20] == 'FE_t0':
		mom_t0_alpha_div.append(species_num_dict[sample[4]])
	elif sample[6][9:11] == 'IS' and sample[6][15:20] == 'FE_t1':
		child_t1_alpha_div.append(species_num_dict[sample[4]])
	elif sample[6][9:11] == 'IS' and sample[6][15:20] == 'FE_t2':
		child_t2_alpha_div.append(species_num_dict[sample[4]])
	elif sample[6][9:11] == 'IS' and sample[6][15:20] == 'FE_t3':
		child_t3_alpha_div.append(species_num_dict[sample[4]])
	elif sample[6][9:11] == 'IS' and sample[6][15:20] == 'FE_t4':
		child_t4_alpha_div.append(species_num_dict[sample[4]])
	elif sample[6][9:11] == 'IS' and sample[6][15:20] == 'FE_t5':
		child_t5_alpha_div.append(species_num_dict[sample[4]])
		
plot_data = [mom_t0_alpha_div, child_t1_alpha_div, child_t2_alpha_div, child_t3_alpha_div, child_t4_alpha_div, child_t5_alpha_div]
plot_names = ["Mothers\nDelivery","1 Day", "3 Days", "1 Week", "1 Month", "4 Months"]

ax[1].boxplot(plot_data)
ax[1].set_xticklabels(plot_names, fontsize=9)
ax[1].set_title("Ferretti")
# plt.legend([bp1["boxes"][i] for i in range(0,4)], plot_names, loc='upper left')

# Yassour
sys.stderr.write("Computing alpha diversity for Yassour\n")

with open(metadata_fpath['yassour'], 'r') as metadata_file:
	metadata = [row.strip().split('\t') for row in metadata_file.readlines()]

shannon_div_dict = {}
species_num_dict = defaultdict(int)
for i in range(len(header)):
	sample = header[i]
	if sample.startswith("SRR728"):
		acc = 0
		for row in data[1:]:
			try:
				rel_ab = float(row[i])
			except ValueError:
				sys.stderr.write('Non-numeric data error\n')
				sys.exit(1)
			if rel_ab != 0:
				acc += (rel_ab * math.log(rel_ab))
				species_num_dict[sample] += 1
		shannon_div_dict[sample] = (acc*-1)

# Unordered lists of all samples' alpha diversities in given cohort
mom_gest_alpha_div = []
mom_birth_alpha_div = []
mom_month3_alpha_div = []

child_birth_alpha_div = []
child_week2_alpha_div = []
child_month1_alpha_div = []
child_month2_alpha_div = []
child_month3_alpha_div = []

for sample in metadata[1:]:
	if sample[7][6:] == 'Mother:Gest':
		mom_gest_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Mother:Birth':
		mom_birth_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Mother:3 months':
		mom_month3_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Child:Birth':
		child_birth_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Child:14 days':
		child_week2_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Child:1 month':
		child_month1_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Child:2 months':
		child_month2_alpha_div.append(species_num_dict[sample[4]])
	elif sample[7][6:] == 'Child:3 months':
		child_month3_alpha_div.append(species_num_dict[sample[4]])
		
plot_data = [mom_gest_alpha_div, mom_birth_alpha_div, mom_month3_alpha_div, child_birth_alpha_div, child_week2_alpha_div, child_month1_alpha_div, child_month2_alpha_div, child_month3_alpha_div]
plot_names = ["Mothers\nGestation","Mothers\nDelivery", "Mothers\n3 Months", "Birth", "2 Weeks", "1 Month", "2 Months", "3 Months"]

ax[2].boxplot(plot_data)
ax[2].set_xticklabels(plot_names, fontsize=9)
ax[2].set_title("Yassour")
# plt.legend([bp1["boxes"][i] for i in range(0,4)], plot_names, loc='upper left')

fig.savefig('%s/alpha_diversity_combined.pdf' % (config.analysis_directory),bbox_inches='tight')
