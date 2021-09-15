import sys, bz2, gzip, os.path
import numpy
from collections import defaultdict
from utils import config, sample_utils as su

# ===========================================================================
# Different possible prevalence cohorts: cohort and directory names
# 
# all			snp_prevalences/all
# infant	snp_prevalences/infant
# hmp			snp_prevalences/hmp
# mother	snp_prevalences/mother
# ===========================================================================

# ===========================================================================
# Standard header to read in argument information
# ===========================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--cohort", help="Cohort to consider for computing prevalence", default="all")
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
cohort = args.cohort

# ===========================================================================

# Open post-processed MIDAS output, annotated_snps.txt
snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)
snp_file =	bz2.BZ2File("%s/annotated_snps.txt.bz2" % snps_dir, 'r')

# Get samples
samples = numpy.array(snp_file.readline().strip().split()[1:])
samples = su.parse_merged_sample_names(samples)

# Filter samples for those in desired cohort
indices = [-1] # For the first column
if cohort == 'hmp':
	hmp_samples = su.get_sample_names('hmp','all')
	for i in range(len(samples)):
		if samples[i] in hmp_samples:
			indices.append(i)
elif cohort == 'infant':
	infant_samples = su.get_sample_names('infant','all')
	for i in range(len(samples)):
		if samples[i] in infant_samples:
			indices.append(i)
elif cohort == 'mother':
	mother_samples = su.get_sample_names('mother', 'all')
	for i in range(len(samples)):
		if samples[i] in mother_samples:
			indices.append(i)
elif cohort == 'premie':
	olm_samples = su.get_sample_names('olm', 'all')
	for i in range(len(samples)):
		if samples[i] in olm_samples:
			indices.append(i)
elif cohort == 'nonpremie':
	infant_samples = su.get_sample_names('infant', 'all')
	olm_samples = su.get_sample_names('olm', 'all')
	for i in range(len(samples)):
		if (samples[i] in infant_samples) and (samples[i] not in olm_samples):
			indices.append(i)
else: # Default is 'all' i.e. leave samples unaltered
	cohort = 'all'

# Get subject-sample map
sample_subject_map = su.parse_sample_subject_map()

record_strs = ["Chromosome, Location, AltFreq, SNPFreq"]

sys.stderr.write("Calculating SNP prevalences...\n")				
num_sites_processed = 0

# Loop through annotated_snps.txt file

for line in snp_file:
		
		num_sites_processed+=1
		
		if num_sites_processed%50000==0:
				sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))		
				if debug:
						break
		
		all_items = line.split()
		
		# Only look at samples in desired prevalence cohort
		items = all_items if cohort == 'all' else [all_items[i+1] for i in indices]
		
		# Load information about site
		info_items = items[0].split("|")
		chromosome = info_items[0]
		location = long(info_items[1])
		gene_name = info_items[2]
		
		# Load alt and depth counts
		
		alts = defaultdict(list) # subject -> list of alt counts
		depths = defaultdict(list) # subject -> list of site depths
		refs = defaultdict(list) # subject -> list of ref counts
		
		for k in range(1, len(items)):
				item = items[k]
				sample = samples[k-1]
				subject = sample_subject_map[sample]
				
				subitems = item.split(",")
				alt, depth = long(subitems[0]), long(subitems[1])
				alts[subject].append(alt)
				depths[subject].append(depth)
				refs[subject].append(depth-alt)
		
		# First calculate fraction of cohort where alternate (non-reference) allele is the major allele in 
		
		population_prevalence = 0
		total_subjects = 0
		
		for subject in depths.keys():
			
			alt_major_samples = 0
			total_samples = 0
			
			# Remember: each subject may have multiple timepoints.
			# As long as site has some coverage in a sample, consider it when
			# computing fraction of samples which have majority alt allele
			
			for i in range(len(depths[subject])):
				depth = depths[subject][i]
				if depth > 0:
					total_samples += 1
					if alts[subject][i] >= refs[subject][i]:
						alt_major_samples += 1
			
			# Define "population prevalence" as sum of proportion of majority-alt
			# samples across all subjects. This is NOT equal to the proportion of
			# all samples which have majority alt allele, because different
			# subjects have different numbers of timepoints. This way prevalence
			# is not biased towards subjects with more timepoints.
			
			if total_samples != 0:
				population_prevalence += float(alt_major_samples)/total_samples
				total_subjects += 1
		
		if total_subjects == 0:
			continue
		else:
			# AKA average proportion of majority-alt samples across subjects
			population_freq = float(population_prevalence)/total_subjects
		
		# Define major allele as having "population frequency" > 50% 
		# If alternate allele is in the majority, it is the major allele
		minors, majors = (refs, alts) if population_freq > 0.5 else (alts, refs)
		
		# Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency
		
		# minor_thresholds = defaultdict(list)
		# for subject in depths.keys():
		# 	for depth in depths[subject]:
		# 		
		# 		alt_threshold = numpy.ceil(depth*0.1)+0.5
		# 		alt_thresholds[subject].append(alt_threshold)
		
		snp_prevalence = 0
		total_subjects = 0
		
		for subject in depths.keys():
			snp_samples = 0
			total_samples = 0
			
			# Remember: each subject may have multiple timepoints.
			# As long as site has some coverage in a sample, consider it when
			# computing fraction of samples which have >=10% minor allele
			
			for i in range(len(depths[subject])):
				depth = depths[subject][i]
				if depth > 0:
					total_samples += 1
					
					# 10% of all reads, rounded up (TODO: why add 0.5?)
					if minors[subject][i] >= numpy.ceil(depth*0.1):
						snp_samples += 1
			
			if total_samples != 0:
				snp_prevalence += float(snp_samples)/total_samples
				total_subjects += 1
		
		snp_freq = float(snp_prevalence)/total_subjects
		
		if (population_prevalence==0) and (snp_prevalence==0):
				continue
		
		# Record nonzero data	 
		record_str = "%s, %d, %g, %g" % (chromosome, location, population_freq, snp_freq) 
		
		record_strs.append(record_str)

snp_file.close()
sys.stderr.write("Done!\n")

# Write to disk!
intermediate_filename_template = config.data_directory + "snp_prevalences/%s/%s.txt.gz"

# Holds panel wide prevalence for each species
os.system('mkdir -p %s/snp_prevalences/%s' % (config.data_directory, cohort))

intermediate_filename = intermediate_filename_template % (cohort, species_name)
output_file = gzip.GzipFile(intermediate_filename,"w")
output_file.write("\n".join(record_strs))
output_file.close
