import sys
from utils import parse_midas_data

if len(sys.argv) > 1:
	species_name=sys.argv[1]
else:
	sys.stderr.write("Usage: pipe_midas_data.py species")

parse_midas_data.pipe_snps(species_name)