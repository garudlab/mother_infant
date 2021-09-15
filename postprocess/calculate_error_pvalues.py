import sys, os, bz2
from utils import config

if len(sys.argv) > 1:
	species_name = sys.argv[1]
else:
	sys.stderr.write("Usage: error_pvalues.py species")

sys.stderr.write("Calculating pvalues for %s...\n" % species_name)

# ===========================================================================
# Creates annotated_snps.txt, which contains alternate and reference allele
# read depths for each sample at each site
# ===========================================================================

snps_dir = "%s/snps/%s/" % (config.data_directory, species_name)
output_filename = "%s/annotated_snps.txt.bz2" % snps_dir

postprocess_dir = '%s/postprocess/' % config.scripts_directory

os.system('python %s/pipe_midas_data.py %s | %s/annotate_pvalue --disabled | bzip2 -c > %s' % (postprocess_dir, species_name, postprocess_dir, output_filename))
