import sys, bz2
import numpy
from collections import defaultdict
from utils import config, sample_utils as su

if len(sys.argv) > 1:
	species_name = sys.argv[1]
else:
	sys.stderr.write("Usage: coverage_distribution.py species")

# Separate into three cohorts: infant, hmp, and adult (including mothers)
prev_cohorts = ['adult', 'infant', 'hmp']

sys.stderr.write("Calculating coverage distribution for %s...\n" % species_name)

# ===========================================================================
# Creates coverage_distribution.txt.bz2 file in species snps directory
# 
# In this file, rows are samples, columns are D,count pairs
# samples are guaranteed to be same order as snps_depth.txt.bz2 file
#
# Also creates a gene_coverage.txt.bz2 file in the species snp directory
# In this file, rows are genes, columns are samples, 
# entries are avg coverage of that gene for that sample 
# ===========================================================================

# These are the default MIDAS parameters. Used to ensure consistency.
# Only sites that pass the prevalence threshold in terms of having at least
# the prevalence_min_coverage can pass, stored in sample_depth_histograms

prevalence_threshold = 0.95
prevalence_min_coverage = 3

# Use all types of sites to include most information
allowed_variant_types = set(["1D","2D","3D","4D"]) 

snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)

# snps_depth.txt: number of reads mapped to genomic site per sample
depth_file = bz2.BZ2File("%s/snps_depth.txt.bz2" % snps_dir, 'r')

# snps_info.txt: metadata for genomic site: mean freq, mean depth, allele props, etc.
info_file = bz2.BZ2File("%s/snps_info.txt.bz2" % snps_dir, 'r')

# Get samples from depth file		
samples = depth_file.readline().split()[1:]

# Get indices of desired prevalence cohort samples
prev_sample_idxs_dict = {}
prev_samples_dict = {}
for prev_cohort in prev_cohorts:
	all_prev_samples = su.get_sample_names(prev_cohort, remove_c=False)
	
	prev_sample_idxs = numpy.array([sample in all_prev_samples for sample in samples])
	prev_sample_idxs_dict[prev_cohort] = prev_sample_idxs
	
	prev_samples = numpy.array(samples)[prev_sample_idxs]
	prev_samples_dict[prev_cohort] = prev_samples

info_file.readline() # remove header

# sample_depth_histograms: (sample -> site -> depth)

# stores only prevalent sites
sample_depth_histograms = {sample: defaultdict(int) for sample in samples}
# stores prevalent sites by prev cohort
sample_depth_histograms_by_prev_cohort = {prev_cohort: {sample: defaultdict(int) for sample in prev_samples_dict[prev_cohort]} for prev_cohort in prev_cohorts}
# stores all sites
full_sample_depth_histograms = {sample: defaultdict(int) for sample in samples}

gene_total_depths = {} # gene -> total depth across all considered sites
gene_total_sites = {} # gene -> number of sites (with known variant type)

num_sites_processed = 0

while True:
		
		# load next lines
		depth_line = depth_file.readline()
		info_line = info_file.readline()
		
		# quit if file has ended
		if depth_line=="":
				break
		
		# parse site info
		info_items = info_line.split('\t')
		variant_type = info_items[5]
				
		# make sure it is a site with known variant type
		if variant_type not in allowed_variant_types:
				continue
		
		gene_name = info_items[6]
		
		depth_items = depth_line.split()
		depths = numpy.array([long(item) for item in depth_items[1:]])
		
		# Manual prevalence filter
		# At least 95% of samples in specified prevalence cohort
		# should have SNP depth >= 3 for this site
		# Compute separately for each prevalence cohort
		
		prevalent_in_some_cohort = False
		
		for prev_cohort in prev_cohorts:
			
			# Restrict samples, depths to those belonging in the prevalence cohort
			prev_sample_idxs = prev_sample_idxs_dict[prev_cohort]
			prev_samples = prev_samples_dict[prev_cohort]
			prev_depths = depths[prev_sample_idxs]
			
			if len(prev_depths) > 0 and (prev_depths>=prevalence_min_coverage).sum()*1.0/len(prev_depths) >= prevalence_threshold: 
					prevalent_in_some_cohort = True
					# Add to genome-wide depth distribution
					for sample, D in zip(prev_samples,prev_depths):
							sample_depth_histograms_by_prev_cohort[prev_cohort][sample][D] += 1
		
		# Alternatively, include the site as long as it has enough SNP depth
		# in 95+% of samples in at least one of the prevalence cohorts,
		# adult/hmp/infant; this list will have no suffix in output file		
		if prevalent_in_some_cohort:
			for sample, D in zip(samples,depths):
				sample_depth_histograms[sample][D] += 1
		
		# Full sample depth histogram (include all sites regardles of prevalence)
		for sample, D in zip(samples,depths):
				full_sample_depth_histograms[sample][D] += 1
		
		# Add to gene-specific avg depth (not subject to prevalence filter)
		if gene_name not in gene_total_depths:
				gene_total_depths[gene_name] = numpy.zeros(len(samples))*1.0
				gene_total_sites[gene_name] = 0.0
				
		gene_total_depths[gene_name] += depths
		gene_total_sites[gene_name] += 1
		
		num_sites_processed += 1
		
		if num_sites_processed%100000==0:
				sys.stderr.write("Processed %dk sites!\n" % (num_sites_processed/1000))

depth_file.close()

# Now write output!

# ===========================================================================
# First write (prevalence-filtered) genome-wide coverage distribution
# Note that sites are included as long as they are prevalent in at least
# one of the three cohorts, adult/hmp/adult
# Note that a site prevalent in hmp may not be prevalent enough in adult...
# I am appending 'any' to not overwrite original coverage_distribution.txt
# 
# Details:
# 	- rows are samples (same as snps_depth.txt samples)
# 	- columns are depths (varies by sample)
# 	- items are number of sites which have that depth, genome-wide
# 
# N.B. we are only considering sites with known variant type, not
# all sites in the genome. Further filtered for having enough coverage
# in at least 95% of all samples.
# ===========================================================================

output_file = bz2.BZ2File("%s/coverage_distribution_any.txt.bz2" % (snps_dir), 'w')
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
		output_file.write("\n")
		output_file.write("\t".join([sample]+["%d,%d" % (D,sample_depth_histograms[sample][D]) for D in sorted(sample_depth_histograms[sample].keys())]))

output_file.close()

# ===========================================================================
# Next write (prevalence-filtered) genome-wide coverage distribution
# where prevalences used for filter are calculated for different cohorts
# ===========================================================================

for prev_cohort in prev_cohorts:
	output_file = bz2.BZ2File("%s/coverage_distribution_prev_%s.txt.bz2" % (snps_dir, prev_cohort), 'w')
	output_file.write("SampleID\tD,n(D) ...")
	prev_samples = prev_samples_dict[prev_cohort]
	for sample in prev_samples:
			output_file.write("\n")
			output_file.write("\t".join([sample]+["%d,%d" % (D,sample_depth_histograms_by_prev_cohort[prev_cohort][sample][D]) for D in sorted(sample_depth_histograms_by_prev_cohort[prev_cohort][sample].keys())]))
	
	output_file.close()

# ===========================================================================
# Write unfiltered genome-wide coverage distribution
# 
# Same form as above, except consider all sites with known variant type.
# ===========================================================================

output_file = bz2.BZ2File("%s/full_coverage_distribution.txt.bz2" % snps_dir, 'w')
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
		output_file.write("\n")
		output_file.write("\t".join([sample]+["%d,%d" % (D, full_sample_depth_histograms[sample][D]) for D in sorted(full_sample_depth_histograms[sample].keys())]))

output_file.close()

# ===========================================================================
# Then write gene-specific coverages
# 
# Details:
# 	- rows are genes
# 	- columns are samples
# 	- items are average depth of all known sites in that gene
# 
# Note: gene_total_sites[gene_name] == 0 should never happen...
# so I got rid of it...
# ===========================================================================

output_file = bz2.BZ2File("%s/gene_coverage.txt.bz2" % snps_dir, 'w')
output_file.write("\t".join(["Gene"]+samples)) # Header line
for gene_name in sorted(gene_total_depths.keys()):
		avg_depths = gene_total_depths[gene_name]/(gene_total_sites[gene_name])
		output_file.write("\n")
		output_file.write("\t".join([gene_name]+["%0.1f" % D for D in avg_depths]))

output_file.close()

# Done!
