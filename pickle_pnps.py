import matplotlib
matplotlib.use('Agg')
from utils import parse_midas_data, sample_utils,  config, sfs_utils, diversity_utils, stats_utils
import sys, os.path, numpy
from math import log10,ceil
import matplotlib.pyplot as plt
from collections import defaultdict
from utils.classes import Interval

type = sys.argv[1] # common, rare, etc.

# Parameters
min_coverage = 20 # orig 20
common_freqrange = [Interval('[0.2, 0.8]')]
rare_freqrange = [Interval('(0, 0.1]'), Interval('[0.9, 1)')]
seg_freqrange = [Interval('(0, 1)')]

if type == 'common':
	freqrange = common_freqrange
elif type == 'rare':
	freqrange = rare_freqrange
elif type == 'seg':
	freqrange = seg_freqrange

# Good species list
good_species_list = parse_midas_data.load_pickled_good_species_list()

# Dictionary: sample -> species -> (within_1D, total_1D, within_4D, total_4D)
within_total_sites_QP = defaultdict(dict)
within_total_sites_nonQP = defaultdict(dict)

# Dictionary: sample -> species -> pN/pS
pNpS_QP = defaultdict(dict)
pNpS_nonQP = defaultdict(dict)

for species in good_species_list:
	
	# Load SNP information for this species
	sys.stderr.write("Loading SFSs for %s...\t" % species)
	samples_4D, sfs_map_4D = parse_midas_data.parse_within_sample_sfs(species, allowed_variant_types=set(['4D'])) # synonymous
	samples_1D, sfs_map_1D = parse_midas_data.parse_within_sample_sfs(species, allowed_variant_types=set(['1D'])) # nonsynonymous
	sys.stderr.write("Done!\n")
	
	# Load genomic coverage distributions
	sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species)
	samples = numpy.array(samples)
	median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(hist) for hist in sample_coverage_histograms])
	sample_median_coverage_map = {samples[i]: median_coverages[i] for i in range(len(samples))}
	
	# Get QP samples (note: low coverage samples are excluded)
	qp_sample_dict = sample_utils.calculate_qp_samples(samples, species)
	samples_QP = qp_sample_dict['qp']
	samples_nonQP = qp_sample_dict['non-qp']
	
	# Only plot samples above a certain median coverage threshold (100)
	desired_samples = samples[(median_coverages >= min_coverage)]
	desired_median_coverages = numpy.array([sample_median_coverage_map[sample] for sample in desired_samples])
	
	if len(desired_samples) <= 0:
		continue
	
	# ====================================
	# Calculate within polymorphism rates
	# ====================================
	
	# Final list of samples used (filtered for nonzero total site counts)
	sample_names = []
	
	# Sites with freq 0-0.05 (rare alleles)
	between_rates_1D = []
	between_rates_4D = []
	
	# Sites with freq = 0, freq > 0.05
	within_rates_1D = []
	within_rates_4D = []
	
	within_sites_1D_array=[]
	total_sites_1D_array=[]
	within_sites_4D_array=[]
	total_sites_4D_array=[]
	
	for sample in desired_samples:
		
		total_sites_1D, within_sites_1D = sfs_utils.calculate_sites_within_freq_range_from_sfs_map(sfs_map_1D[sample], freqrange)
		total_sites_4D, within_sites_4D = sfs_utils.calculate_sites_within_freq_range_from_sfs_map(sfs_map_4D[sample], freqrange)
		
		# Skip if zero of either syn. or nonsyn. total sites
		if total_sites_1D <= 0 or total_sites_4D <= 0:
			continue
		
		# Fraction of all nonsynonymous sites with minor allele frequency > 0.05
		pN = (within_sites_1D*1.0 + 1.0)/(total_sites_1D + 1.0)		
		# Fraction of all synonymous sites with minor allele frequency > 0.05
		pS = (within_sites_4D*1.0 + 1.0)/(total_sites_4D + 1.0)
		
		# Store within and total sites, pN/pS for each sample-species pair
		if sample in samples_QP:
			within_total_sites_QP[sample][species] = (within_sites_1D, total_sites_1D, within_sites_4D, total_sites_4D)
			pNpS_QP[sample][species] = pN/pS
		
		elif sample in samples_nonQP:
			within_total_sites_nonQP[sample][species] = (within_sites_1D, total_sites_1D, within_sites_4D, total_sites_4D)
			pNpS_nonQP[sample][species] = pN/pS

# Pickle!!
import pickle

pdir = "%s/pickles" % config.data_directory

pickle.dump(within_total_sites_nonQP, open("%s/within_total_sites_%s_nonQP_cov20_rare10pct.pkl" % (pdir, type), 'wb'))
pickle.dump(within_total_sites_QP, open("%s/within_total_sites_%s_QP_cov20_rare10pct.pkl" % (pdir, type), 'wb'))
pickle.dump(pNpS_nonQP, open("%s/pNpS_%s_nonQP.pkl" % (pdir, type), 'wb'))
pickle.dump(pNpS_QP, open("%s/pNpS_%s_QP.pkl" % (pdir, type), 'wb'))
