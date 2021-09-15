# Script to pickle very specific information (allele identity and counts)
# for small set of given QP pairs at given sites

from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils, core_gene_utils, gene_diversity_utils
import numpy as np
from numpy.random import choice, random as np_random, randint
import random
from collections import defaultdict
import pickle
import bz2
import numpy
import os, sys

# Desired samples and sites
# sites are given as gene_id, contig, location tuples
# inconsistent but bleh

desired_host_species_sites = [(['M0806-C'], \
  'Bacteroides_vulgatus_57955', \
  [('435590.9.peg.1499', 'NC_009614', 1915720L), \
   ('435590.9.peg.1499', 'NC_009614', 1915720L)]), \
 (['3-I'], \
  'Bifidobacterium_catenulatum_58257', \
  [('566552.4.peg.96', 'NZ_ABXY01000001', 124565L), \
   ('566552.4.peg.96', 'NZ_ABXY01000001', 124565L)]), \
 (['67-I'], \
  'Bacteroides_fragilis_54507', \
  [('1339327.3.peg.4421', 'JGDJ01000264', 9480L), \
   ('1339327.3.peg.4421', 'JGDJ01000264', 9480L)]), \
 (['M0901-C'], \
  'Blautia_wexlerae_56130', \
  [('1121115.4.peg.32', 'AXVN01000001', 35625L), \
   ('1121115.4.peg.32', 'AXVN01000001', 35625L)]), \
 (['1-I'], \
  'Bacteroides_fragilis_54507', \
  [('1339327.3.peg.2283', 'JGDJ01000171', 5687L), \
   ('1339327.3.peg.2283', 'JGDJ01000171', 5687L)])]

# desired_host_species_sites = [(['C02143-I', 'C02143-M'], 'Bifidobacterium_bifidum_55065', [('500634.3.peg.1861','AWSW01000054',37945), ('500634.3.peg.945','AWSW01000030',35925), ('500634.3.peg.1636','AWSW01000046',21960), ('500634.3.peg.875','AWSW01000027',7916), ('500634.3.peg.952','AWSW01000030',45351), ('500634.3.peg.1619','AWSW01000046',4226)]), \
# (['59-I', '59-M'], 'Bifidobacterium_adolescentis_56815', [('592977.3.peg.642','JGZQ01000005',14361),('592977.3.peg.860','JGZQ01000006',69020),('592977.3.peg.1216','JGZQ01000008',284119),('592977.3.peg.1129','JGZQ01000008',186577), ('592977.3.peg.39','JGZQ01000001',53358), ('592977.3.peg.1732','JGZQ01000009',58203), ('592977.3.peg.1705','JGZQ01000009',29026)])]

# desired_samples = ['ERR3405741', 'ERR3405661', 'ERR3406235']

# ===========================================================================
# Loads allele counts for specific samples at specific sites
# where sites are provided as (contig, location, gene_id) tuples
# TODO: move to parse_midas_data later
# ===========================================================================

def parse_snps_specify_sites_details(species_name, desired_samples=[], desired_sites=[], prev_cohort='all'):
	
	# Alternate version without gene names
	desired_sites_no_gene = [(contig, location) for contig, location, gene_id in desired_sites]
	
	# SNPs directory
	snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)
	
	# Load population freqs (for polarization purposes)
	population_freqs = snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=False)
	
	# Open post-processed MIDAS output
	snp_file = bz2.BZ2File("%s/annotated_snps.txt.bz2" % snps_dir, 'r')
	
	# =====================================================================
	# Process allele information
	# =====================================================================	
	
	# sample -> site -> (ref allele, alt allele)
	sample_site_allele_dict = defaultdict(dict)
	
	# Load snps_alt_allele.txt
	snp_alleles_file = bz2.BZ2File("%s/snps_alt_allele.txt.bz2" % snps_dir, 'r')
	items = snp_alleles_file.readline().strip().split()[1:]
	samples = list(su.parse_merged_sample_names(items))
	
	desired_sample_idxs = []
	for sample in desired_samples:
		if sample in samples:
			desired_sample_idxs.append(samples.index(sample))
	
	desired_sample_idxs = numpy.array(sorted(desired_sample_idxs))
	
	for line in snp_alleles_file:		
		items = line.split()
		
		# Load information about site
		info_items = items[0].split("|")
		contig = info_items[0]
		location = long(info_items[1])
		ref_allele = info_items[2]
		
		if (contig, location) not in desired_sites_no_gene:
			continue
		
		for idx in desired_sample_idxs:
			alt_allele = items[1+idx]
			sample = samples[idx]
			sample_site_allele_dict[sample][(contig, location)] = (ref_allele, alt_allele)
	
	snp_alleles_file.close()
	
	# =====================================================================
	# Process annotated_snps information
	# =====================================================================
	
	# Open post-processed MIDAS output	
	snp_file = bz2.BZ2File("%s/annotated_snps.txt.bz2" % snps_dir, 'r')
	
	# Get lists of desired sample idxs
	items = snp_file.readline().strip().split()[1:]
	samples = list(su.parse_merged_sample_names(items))
	
	desired_sample_idxs = []
	for sample in desired_samples:
		if sample in samples:
			desired_sample_idxs.append(samples.index(sample))
	
	desired_sample_idxs = numpy.array(sorted(desired_sample_idxs))
	
	# Map: sample -> site (contig, location, gene_id) -> allele count
	allele_counts_map = defaultdict(dict)
	
	# Map: site (contig, location, gene_id) -> variant type
	variant_type_map = {}
	
	num_sites_processed = 0
	
	# Loop over sites in annotated_snps.txt file
	for line in snp_file:
		if num_sites_processed>0 and num_sites_processed%10000==0:
			sys.stderr.write("%d0k sites processed...\n" % (num_sites_processed/10000))
		num_sites_processed += 1
		
		items = line.split()
		
		# Load information about site
		info_items = items[0].split("|")
		contig = info_items[0]
		location = long(info_items[1])
		gene_name = info_items[2]
		variant_type = info_items[3]
		polarization = 'R' # note R was assigned indiscriminately
		pvalue = float(info_items[5])
		
		# Only look at sites of interest
		if (contig, location, gene_name) not in desired_sites:
			continue
		
		# Store variant type
		variant_type_map[(contig, location, gene_name)] = variant_type
		
		# Store alt and depth counts at this site for all desired samples
		for idx in desired_sample_idxs:
			alt, depth = [float(num) for num in items[1+idx].split(",")]
			sample = samples[idx]
			allele_counts_map[sample][(contig, location, gene_name)] = (alt, depth)
	
	snp_file.close()
	
	return allele_counts_map, sample_site_allele_dict, variant_type_map

# Load a few things
subject_sample_map = su.parse_subject_sample_map()

# Set up pickle directory

ddir = config.data_directory
pdir = "%s/pickles/reversion_examples/" % (ddir)
os.system('mkdir -p %s' % pdir)

# Store these two dicts for each host-species pair
for subjects, species, sites in desired_host_species_sites:
	
	# Get all samples within the host
	desired_samples = []
	for subject in subjects:
		desired_samples += subject_sample_map[subject].keys()
	
	# Reformat sites
	desired_sites = []
	for gene_id, contig, location in sites:
		desired_sites.append((contig, location, gene_id))
	
	# Obtain desired dicts
	allele_counts_map, sample_site_allele_dict, variant_type_map = parse_snps_specify_sites_details(species, desired_samples, desired_sites=desired_sites)
	
	# Pickle dicts
	sys.stderr.write("Pickling...\n")
	pickle.dump((allele_counts_map, sample_site_allele_dict, variant_type_map), open('%s/allele_info_%s_%s.pkl' % (pdir, ('_').join(subjects), species), 'wb'))

sys.stderr.write("Done!\n")
