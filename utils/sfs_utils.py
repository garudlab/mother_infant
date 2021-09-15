##############################################################################
#
# Utility methods for working with (primarily within-sample) SFSs
#
##############################################################################
import numpy
import config
from math import log10
from classes import Interval

def calculate_binned_sfs_from_sfs_map(sfs_map, bins=[], folding='minor'):

    alts = []
    depths = []
    counts = []
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        
        alts.append(A)
        depths.append(D)
        counts.append(n)
    
    alts = numpy.array(alts)
    depths = numpy.array(depths)
    counts = numpy.array(counts)
    weights = counts*1.0/counts.sum()
    
    freqs = alts*1.0/depths
    minor_freqs = numpy.fmin(freqs,1-freqs)
    
    # calculate median depth (or rough approximation)
    sorted_depths, sorted_weights = (numpy.array(x) for x in zip(*sorted(zip(depths, weights))))
    CDF = numpy.cumsum(sorted_weights)
    Dbar = sorted_depths[CDF>0.5][0]
    
    if len(bins)==0:
        # use this to set up bins
        bins = (numpy.arange(0,Dbar+2)-0.5)/Dbar
        fs = numpy.arange(0,Dbar+1)*1.0/Dbar
        
    else:
        bins = numpy.array(bins)
        fs = bins[1:]
    
    pfs = numpy.zeros_like(fs)
    bin_idxs = numpy.digitize(minor_freqs, bins=bins)
    
    for bin_idx, weight in zip(bin_idxs, weights):
        pfs[bin_idx-1] += weight
    
    # should already be normalized, but just to make sure...
    pfs /= pfs.sum()
    
    if folding=='major':
        pfs = pfs[::-1]
        fs = (1.0-fs)[::-1]
    
    return fs,pfs
    
def calculate_binned_depth_distribution_from_sfs_map(sfs_map,bins=[],num_bins=30):

    alts = []
    depths = []
    counts = []
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        
        alts.append(A)
        depths.append(D)
        counts.append(n)
    
    alts = numpy.array(alts)
    depths = numpy.array(depths)
    counts = numpy.array(counts) # read counts for all sites with (alt, depth)
    weights = counts*1.0/counts.sum() # read counts for each frequency bin normalized by total read counts
    
    # calculate median depth (or rough approximation)
		
		# sort weights by depth
    sorted_depths, sorted_weights = (numpy.array(x) for x in zip(*sorted(zip(depths, weights))))
		
		# get median depth
    CDF = numpy.cumsum(sorted_weights)
    Dbar = sorted_depths[CDF>0.5][0]
    
    if len(bins)==0:
        # use this to set up bins
				# range from median depth/8 to depth*8, log space
        bins = numpy.logspace(log10(Dbar/8),log10(Dbar*8),num_bins)
				# Ds are depth bins
        Ds = bins[0:-1]
           
    else:
        bins = numpy.array(bins, copy=True)
        Ds = bins[0:-1]
    
    bins[0] = 0
    bins[-1] = 1e09
    
    pDs = numpy.zeros_like(Ds)
		# For each depth, assign which bin it falls in (via bins index)
    bin_idxs = numpy.digitize(depths, bins=bins)
    
    for bin_idx, weight in zip(bin_idxs, weights):
				# pDs are weights for depth bins
        pDs[bin_idx-1] += weight
    
    # should already be normalized, but just to make sure...
    pDs /= pDs.sum()
    
    return bins, Ds, pDs

def calculate_depth_distribution_from_sfs_map(sfs_map):
    
    depth_map = {}
    for key in sfs_map.keys():
        D,A = key
        count = sfs_map[key][0]
        
        if D not in depth_map:
            depth_map[D] = 0
        
        depth_map[D] += count
        
    depths = numpy.array(sorted(depth_map.keys()))
    counts = numpy.array([depth_map[d] for d in depths])
    
    return depths, counts

def calculate_polymorphism_rates_from_sfs_map(sfs_map,lower_threshold=0.2,upper_threshold=0.8):
    
    total_sites = 0 # Total number of sites
    within_sites = 0 # Number of sites with frequency 0.2-0.8 (not inclusive)
    between_sites = 0 # Number of sites with frequency <= 0.2, >= 0.8
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        reverse_n = sfs_map[key][1]
        
        f = A*1.0/D
        
        total_sites += n
        
        if ((f>lower_threshold) and (f<upper_threshold)):
            # an intermediate frequency site
            within_sites += n
        else:    
            if f>0.5:
                between_sites += (n-reverse_n) #NRG what does this do?
            else:
                between_sites += reverse_n
        
        
    return within_sites, between_sites, total_sites

# =============================================================
# number of polymorphisms among synonymous sites
# (number sites with minor allele frequency >0.05) / (total number synonymous sites)
# Rare SFS: (0, 0.05], [0.95, 1)
# Common SFS: [0.2, 0.8]
# =============================================================

# =============================================================
# Computes number of sites, in given SFS, that have frequency
# within given list of frequency intervals/ranges
# 
# Note that Interval class is used to express ranges
# =============================================================

def calculate_sites_within_freq_range_from_sfs_map(sfs_map, freq_ranges = [Interval('[0.2, 0.8]')], min_alt = 2):
	
	total_sites = 0
	within_sites = 0
	
	for D, A in sfs_map:
		n, reverse_n = sfs_map[(D, A)]
		if A < min_alt:
			f = 0 # If not enough alt reads, assume frequency 0
		else:
			f = A*1.0/D # Frequency of alternate allele
		
		total_sites += n
		
		for range in freq_ranges:
			# As long as one of the ranges contains the frequency
			if range.contains(f):
				within_sites += n
				continue
	
	return total_sites, within_sites

def calculate_singleton_rates_from_sfs_map(sfs_map, lower_threshold=0, upper_threshold=0.2):
  
	total_sites = 0 # Total number of sites in SFS
	within_sites = 0 # Frequency < lower_threshold, > upper_threshold
	between_sites = 0 # Frequency within lower_threhsold-upper_threshold
	
	for D, A in sfs_map.keys():
		n, reverse_n = sfs_map[(D, A)]
		f = A*1.0/D # Frequency of alternate allele
		
		total_sites += n
		
		if ((f>0 and f<lower_threshold) or (f>upper_threshold and f<1)):
			# an intermediate frequency site
			within_sites += n
		else:    
			if f>0.5:
					between_sites += (n-reverse_n) #NRG what does this do?
			else:
					between_sites += reverse_n
	
	return within_sites, between_sites, total_sites
