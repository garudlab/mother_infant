from collections import defaultdict

coverage_threshold = 3.0

Olm_campaign = 'NIH2'

Olm_full_directory = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Olm_2019/%s' % Olm_campaign
metadata_directory = '/u/project/ngarud/daisyche/scripts/metadata'

sample_ids_file = open('%s/Olm_%s_run_accessions.txt' % (Olm_full_directory, Olm_campaign), 'r')
metadata_file = open('%s/olm_accessions_babyIds.txt' % metadata_directory, 'r')
merge_file = open('%s/Olm_%s_merge_replicates.txt' % (Olm_full_directory, Olm_campaign), 'r')

# Get all relevant sample IDs for this campaign
all_run_accs = []
for line in sample_ids_file:
	sample = line.strip()
	all_run_accs.append(sample)

# Get merge information
acc_comb_merge_dict = {}
comb_acc_merge_dict = {}
for line in merge_file:
	acc1, acc2, comb = line.split()
	acc_comb_merge_dict[acc1] = comb
	acc_comb_merge_dict[acc2] = comb
	comb_acc_merge_dict[comb] = [acc1, acc2]

# Maps subject ID to dictionary of timepoint-run accession
subject_sample_time_map = defaultdict(dict)

metadata_file.readline() # Ignore first line

for line in metadata_file:
	try:
		public_code, subject_id, days, campaign, _, _, _, _, _, run_acc = line.strip().split()
		timept = days
		if run_acc in all_run_accs:
			if run_acc in acc_comb_merge_dict:
				comb_acc = acc_comb_merge_dict[run_acc]
				subject_sample_time_map[subject_id][timept] = comb_acc
			else:
				subject_sample_time_map[subject_id][timept] = run_acc
	except:
		pass # No run acc available

bad_samples = [] # Run accessions for samples lacking species_profile.txt

samples_in_subject_with_merge = []

# Figure out which samples are in subjects that have a merged sample
for subject_id in subject_sample_time_map:
	for timept in subject_sample_time_map[subject_id]:
		run_acc = subject_sample_time_map[subject_id][timept]
		if run_acc in comb_acc_merge_dict:
			samples_in_subject_with_merge += subject_sample_time_map[subject_id].values()

# Get species union set for each subject
for subject_id in subject_sample_time_map:
	
	species_union = set()
	
	for timept in subject_sample_time_map[subject_id]:
		run_acc = subject_sample_time_map[subject_id][timept]
		
		# Get species list for one sample
		species_list = set()
		
		try:
			midas_species_file = open('%s/midas_output/%s/species/species_profile.txt' % (Olm_full_directory, run_acc), 'r')
		except:
			bad_samples.append(run_acc)
			break
		
		midas_species_file.readline() # Ignore first line
		for line in midas_species_file:
			species_id, _, coverage, _ = line.strip().split()
			if float(coverage) >= coverage_threshold: 
				species_list.add(species_id)
		
		# Add to the species union set
		species_union = species_union.union(species_list)
	
	# Now species union should have species from all timepoints for this subject
	
	for timept in subject_sample_time_map[subject_id]:
		run_acc = subject_sample_time_map[subject_id][timept]
		
		# Species union output file
		species_union_file = open('%s/species_unions/%s_species_union.txt' % (Olm_full_directory, run_acc), 'w') 
		species_union_file.write('species_id\n') # header
		for species in species_union:
				species_union_file.write(species + '\n')
