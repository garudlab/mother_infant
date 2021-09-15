# Question: are sweeping alleles in infants present in mom?
# Idea: store data on allele frequency in mom for each sweeping allele in infant

from utils import sample_utils as su, parse_midas_data, substitution_rates_utils, config, temporal_changes_utils, snps_utils, core_gene_utils, gene_diversity_utils
import numpy as np
from numpy.random import choice, random as np_random, randint
import random
from collections import defaultdict
import pickle
import bz2
import numpy
import os, sys

# ===========================================================================
# Loads allele counts for specific samples at specific sites
# where sites are provided as (contig, location, gene_id) tuples
# TODO: move to parse_midas_data later
# ===========================================================================

def parse_snps_specify_sites(species_name, desired_samples=[], desired_sites=[], prev_cohort='all'):
	
	# Load population freqs (for polarization purposes)
	population_freqs = snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=False)
	
	# Open post-processed MIDAS output
	snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)
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
		polarization = 'R' # note R was assigned indiscriminately
		pvalue = float(info_items[5])
		
		# Only look at sites of interest
		if (contig, location, gene_name) not in desired_sites:
			continue
		
		# Load alt and depth counts at this site for all desired samples
		alts, depths = [], []
		
		for idx in desired_sample_idxs:
			alt, depth = [float(num) for num in items[1+idx].split(",")]
			alts.append(alt); depths.append(depth)
		
		alts = numpy.array(alts)
		depths = numpy.array(depths)
		
		# Obtain population frequency of alt allele
		# Recall: this is average proportion of majority-alt samples across subjects
		if (contig, location) in population_freqs:
			population_freq = population_freqs[(contig, location)]
		else: # alt population prevalence is (probably? TODO) 0
			population_freq = 0
		
		# Polarize SFS according to population freq
		if population_freq > 0.5: # This means alt allele is the major allele
			alts = depths - alts
			polarization = 'A'
		
		# For sites from temporal changes, note that we should assume
		# site is "passed" i.e. alt aleles make up >5% of all alleles
		# in at least one sample and pvalue < 0.05
		
		# Store allele counts only if the site is interesting
		for alt, depth, sample in zip(alts, depths, desired_samples):
			allele_counts_map[sample][(contig, location, gene_name)] = (alt, depth)
	
	snp_file.close()
	
	return allele_counts_map

# Parameters
species = sys.argv[1]
sweep_type = 'full' # assume full for now
pp_prev_cohort = 'all'
min_coverage = 0
thresholds = {'full': (0.2, 0.8), 'partial': (0.35, 0.65)}
lower_threshold, upper_threshold = thresholds[sweep_type]

# Sample-subject-order maps
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_subject_map = su.parse_sample_subject_map()
sample_cohort_map = su.parse_sample_cohort_map()
same_mi_pair_dict = su.get_same_mi_pair_dict(subject_sample_map)
sys.stderr.write("Done!\n")

# Use output from pickle_everything.py to get list of sweeping alleles in infants
ddir = config.data_directory
pdir = "%s/pickles/cov%i_prev_%s/nonconsecutive/" % (ddir, min_coverage, pp_prev_cohort)
snp_changes = pickle.load(open('%s/big_snp_changes_%s.pkl' % (pdir, sweep_type), 'rb'))
'''
# Loop over all modification SNP changes
for species in snp_changes:
	
	# Form list of sites and samples of interest
	desired_sites = set()
	desired_samples = set()
	
	for s1, s2 in snp_changes[species]:
		
		val = snp_changes[species][(s1, s2)]
		
		# Interested in all samples within the mother-infant pair
		# Note that mothers and infants differ by 2-letter suffix of subject
		subject = sample_subject_map[s1]		
		for sample in subject_sample_map[subject]:
			desired_samples.add(sample)
		
		if subject in same_mi_pair_dict:
			other_subject = same_mi_pair_dict[subject]
			for sample in subject_sample_map[other_subject]:
				desired_samples.add(sample)
		
		if type(val) == type([]): # Modification event
			for snp_change in val:
				gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
				desired_sites.add((contig, position, gene_name))
	
	# Load allele_counts_map
	allele_counts_map = parse_snps_specify_sites(species, desired_samples, desired_sites=list(desired_sites))
'''
# Form list of sites and samples of interest
desired_sites = set()
desired_samples = set()

for s1, s2 in snp_changes[species]:
	
	val = snp_changes[species][(s1, s2)]
	
	# Only look at mothers and infants excluding olm
	cohort = sample_cohort_map[s1]
	if cohort not in ['backhed', 'ferretti', 'yassour', 'shao']:
		continue
	
	# Interested in all samples within the mother-infant pair
	# Note that mothers and infants differ by 2-letter suffix of subject
	subject = sample_subject_map[s1]		
	for sample in subject_sample_map[subject]:
		desired_samples.add(sample)
	
	if subject in same_mi_pair_dict:
		other_subject = same_mi_pair_dict[subject]
		for sample in subject_sample_map[other_subject]:
			desired_samples.add(sample)
	
	if type(val) == type([]): # Modification event
		for snp_change in val:
			gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
			desired_sites.add((contig, position, gene_name))

# Load allele_counts_map
allele_counts_map = parse_snps_specify_sites(species, desired_samples, desired_sites=list(desired_sites))

# Pickle time
sys.stderr.write("Pickling...\n")

ddir = config.data_directory
pdir = "%s/pickles/seeding/" % (ddir)
os.system('mkdir -p %s' % pdir)

pickle.dump(allele_counts_map, open('%s/allele_counts_map_%s.pkl' % (pdir, species), 'wb'))

sys.stderr.write("Done!\n")

