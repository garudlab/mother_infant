from collections import defaultdict

coverage_threshold = 3.0

Shao_directory = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/Shao_2019'
metadata_directory = '/u/project/ngarud/daisyche/mother_infant/scripts/metadata'

sample_ids_file = open('%s/PRJEB32631_run_accession_only.txt' % Shao_directory, 'r')
metadata_file = open('%s/Shao_ids.txt' % metadata_directory, 'r')

# Maps subject ID to dictionary of timepoint-run accession
subject_sample_time_map = defaultdict(dict)

metadata_file.readline() # Ignore first line

for line in metadata_file:
	run_acc, _, subject_id, days, infancy_months = line.strip().split()
	timept = "I" + infancy_months if days == "Infancy" else days
	subject_sample_time_map[subject_id][timept] = run_acc

bad_samples = [] # Run accessions for samples lacking species_profile.txt

# Get species union set for each subject
for subject_id in subject_sample_time_map:
	
	species_union = set()
	
	for timept in subject_sample_time_map[subject_id]:
		run_acc = subject_sample_time_map[subject_id][timept]
		
		# Get species list for one sample
		species_list = set()
		
		try:
			midas_species_file = open('%s/midas_output/%s/species/species_profile.txt' % (Shao_directory, run_acc), 'r')
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
		species_union_file = open('%s/species_unions/%s_species_union.txt' % (Shao_directory, run_acc), 'w') 
		species_union_file.write('species_id\n') # header
		for species in species_union:
				species_union_file.write(species + '\n')
