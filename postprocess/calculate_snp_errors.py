# Within-snp gene changes

from utils import config, parse_midas_data, sample_utils, diversity_utils, gene_diversity_utils, sfs_utils, stats_utils
import sys
import numpy

species_name = "Bacteroides_vulgatus_57955"

# =================================================
# Standard header to read in argument information
# =================================================
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""

# =================================================

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
clade_divergence_threshold = 1e-02
modification_divergence_threshold = 1e-03
min_change = 0.8
include_high_copynum = False
#include_high_copynum = True

sample_size = 3e06

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_order_map = sample_utils.parse_sample_order_map()
sys.stderr.write("Done!\n")

temporal_samples = diversity_utils.calculate_temporal_samples(species_name, min_coverage=config.min_median_coverage)
haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name, debug=debug))


from utils import sfs_utils
sys.stderr.write("Loading SFSs for %s...\t" % species_name)
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
sys.stderr.write("Done!\n")

same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, temporal_samples)

frequency_bins = numpy.linspace(0,1,21)

#fs = numpy.array([0.1,0.2,0.3,0.4,0.5,0.5,0.6,0.7,0.8,0.9])

dfs = numpy.array([0.6,0.7,0.8,0.9,0.98])

perrs = []

for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
       
    i = same_subject_idxs[0][sample_pair_idx]
    j = same_subject_idxs[1][sample_pair_idx]

    sample_i = temporal_samples[i]
    sample_j = temporal_samples[j]
    
    if sample_i in haploid_samples:
        ploidy_i = 'haploid'
    else:
        ploidy_i = 'diploid'
        
    if sample_j in haploid_samples:
        ploidy_j = 'haploid'
    else:
        ploidy_j = 'diploid'
    
    # Calculate SFS distributions    
    dummy_fs, pfs_i = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_i],bins=frequency_bins)
    dummy_fs, pfs_j = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_j],bins=frequency_bins)
    
    fs = frequency_bins[1:]-(frequency_bins[1]-frequency_bins[0])/2.0

    pfs = (pfs_i+pfs_j)/2.0 # Average out SFS distribution    
    pfs = (pfs+pfs[::-1])/2 # Fold (make symmetric)
    
    # Calculate depth distributions
    bins1, D1s, pD1s = sfs_utils.calculate_binned_depth_distribution_from_sfs_map(sfs_map[sample_i])
    bins2, D2s, pD2s = sfs_utils.calculate_binned_depth_distribution_from_sfs_map(sfs_map[sample_j])
    		
    fs = fs[pfs>0] 
    pfs = pfs[pfs>0]
    
		# Restrict to bins with nonzero weight
    D1s = D1s[pD1s>0]
    pD1s = pD1s[pD1s>0]
    D2s = D2s[pD2s>0]
    pD2s = pD2s[pD2s>0]
    
    from scipy.stats import binom
    
    perrs.append({df:0 for df in dfs})
    for D1,pD1 in zip(D1s,pD1s):
        for D2,pD2 in zip(D2s,pD2s):
            for f,pf in zip(fs,pfs):
                for df in dfs:
									k1 = D1*(1-df)/2 # Minor allele count
									n1 = D1 # Observed depth [bin]
									p1 = f # True minor allele frequency
									
									k2 = D2*(1-df)/2
									n2 = D2
									p2 = 1-f
									
                  perrs[-1][df] += binom.cdf(k1, n1, p1)*binom.cdf(k2, n2, p2)*pD1*pD2*pf
    
    if True: #ploidy_i=='haploid' and ploidy_j=='haploid':
        print sample_i, D1s[len(D1s)/2], ploidy_i, sample_j, ploidy_j, D2s[len(D2s)/2], perrs[-1][0.6]*3e06    
        
        
print "Done!"
        
    
        
