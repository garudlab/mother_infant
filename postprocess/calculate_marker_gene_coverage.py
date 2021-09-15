import sys, bz2
import numpy
from utils import config

if len(sys.argv) > 1:
	species_name = sys.argv[1]
else:
	sys.stderr.write("Usage: marker_gene_coverage.py species")

sys.stderr.write("Calculating marker gene coverage for %s...\n" % species_name)

# ===========================================================================
# Creates species-specific marker_coverage.txt.bz2 in species snps directory
# ===========================================================================

snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)
species_dir = "%s/species" % config.data_directory

# Get list of samples in depth file
depth_file = bz2.BZ2File("%s/snps_depth.txt.bz2" % snps_dir, 'r')
output_samples = depth_file.readline().split()[1:]
depth_file.close()

# coverage.txt: average read depth of 15 marker genes per species
coverage_file = bz2.BZ2File("%s/coverage.txt.bz2" % species_dir, 'r')

# Get list of samples in coverage file (should be full list)
samples = coverage_file.readline().strip().split()[1:]

# Get the indices of depth samples within the list of coverage file samples
output_idxs = numpy.array([samples.index(sample) for sample in output_samples])

# marker_coverage.txt: 
output_coverage_file = bz2.BZ2File("%s/marker_coverage.txt.bz2" % snps_dir, 'w')
	
# write header lines for output file
output_coverage_file.write('\t'.join(['species_id'] + output_samples) + '\n')

total_depths = numpy.zeros(len(output_samples))

for line in coverage_file:
	
	items = line.split()
	current_species = items[0]
	depths = numpy.array([float(d) for d in items[1:]])[output_idxs]
	
	if current_species == species_name:
		output_coverage_file.write("\t".join([current_species]+["%g" % d for d in depths]))
	
	# Total depths reports sum of average marker depth across all species
	total_depths += depths

output_coverage_file.write('\n')
output_coverage_file.write('\t'.join(["Total"] + ["%g" % d for d in total_depths]))

coverage_file.close()
output_coverage_file.close()