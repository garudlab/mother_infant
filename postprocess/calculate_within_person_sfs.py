import sys, bz2
import numpy
from utils import config, parse_midas_data, snps_utils

# ====================================================================================
# Standard header to read in argument information
# ====================================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--cohort", help="Cohort to consider for computing prevalence", default="all")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
args = parser.parse_args()

species_name = args.species_name
cohort = args.cohort
debug = args.debug

# ====================================================================================

# Only look at core genes
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name, prev_cohort='all')
sys.stderr.write("Done! %d core genes\n" % len(core_genes))
allowed_genes = core_genes

# Load population frequencies
sys.stderr.write("Loading population freqs...\n")
population_freqs = snps_utils.parse_population_freqs(cohort, species_name)
sys.stderr.write("Done! %d SNVs\n" % len(population_freqs))

# Allowed variant types
allowed_variant_types = set(['1D','2D','3D','4D'])
allowed_variant_type_list = ['1D','2D','3D','4D']

# Open post-processed MIDAS snps output
snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)
snp_file =	bz2.BZ2File("%s/annotated_snps.txt.bz2" % snps_dir, 'r')

# Remove header and get samples
samples = numpy.array(snp_file.readline().strip().split()[1:])

# Stores sites by variant type for each sample
site_map = [{vartype: {} for vartype in allowed_variant_types} for _ in samples]

sys.stderr.write("Calculating within-person SFSs...\n")				 
num_sites_processed = 0

# Examine each site (in all samples), one by one

for line in snp_file:
	
	items = line.split()
	
	# Load information about site
	
	info_items = items[0].split("|")
	chromosome = info_items[0]
	location = long(info_items[1])
	gene_name = info_items[2]
	variant_type = info_items[3]
	
	if len(info_items) > 5: # for backwards compatability
		polarization = info_items[4]
		pvalue = float(info_items[5])
	else: 
		polarization="?"
		pvalue = float(info_items[4])
	
	# Ignore site if variant type not known, or if it's in a bad gene
	
	if variant_type not in allowed_variant_types:
		continue
	
	if gene_name not in allowed_genes:
		continue
	
	# Load alt and depth counts
	
	alts = []
	depths = []
	
	for item in items[1:]:
		subitems = item.split(",")
		alts.append(long(subitems[0]))
		depths.append(long(subitems[1]))
	
	alts = numpy.array(alts)
	depths = numpy.array(depths)
	refs = depths - alts
	
	# population_freq returns the fraction of people for which the alt is the major allele.
	# This is a very important quantity being computed! It is later used for identifying CPS samples.
	if (chromosome, location) in population_freqs:
			population_freq = population_freqs[(chromosome, location)]
	else: # This is because sites are skipped when pop. prevalence (or snp_prevalence) is 0
			population_freq = 0
	
	# Polarize SFS according to population frequency of alternate allele
	# which we now define as population frequency of minor allele (could be ref)
	
	minors, majors = (refs, alts) if population_freq > 0.5 else (alts, refs)
	population_freq = (1 - population_freq) if population_freq > 0.5 else population_freq
	
	for i in range(len(minors)): # i indexes the sample
		
		site = (depths[i], minors[i]) # indicates proportion of minor alleles
		
		if site not in site_map[i][variant_type]:
			site_map[i][variant_type][site] = [0,0.0]
		
		# For each (depth,minor) in each sample, report number of sites that have this
		# read distribution, and report minor allele population frequency
		
		site_map[i][variant_type][site][0] += 1
		site_map[i][variant_type][site][1] += population_freq # weight of polarization reversals
	
	num_sites_processed+=1
	
	if num_sites_processed % 50000 == 0:
			sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))		
			if debug:
					break

snp_file.close()
sys.stderr.write("Done!\n")

# Write to disk!
sys.stderr.write("Writing output...\n")

# First write (filtered) genome-wide coverage distribution
output_file = bz2.BZ2File("%s/within_sample_sfs_prev_%s.txt.bz2" % (snps_dir, cohort), 'w')
output_file.write("\t".join(["SampleID", "variant_type", "D,A,count,reverse_count", "..."]))

for sample_idx in range(len(samples)):
	sample = samples[sample_idx]
	for variant_type in allowed_variant_type_list:
		output_file.write("\n")
		output_file.write("\t".join([sample, variant_type]+["%d,%d,%d,%g" % (site[0],site[1],site_map[sample_idx][variant_type][site][0],site_map[sample_idx][variant_type][site][1]) for site in sorted(site_map[sample_idx][variant_type].keys())]))

output_file.close()
sys.stderr.write("Done!\n")
