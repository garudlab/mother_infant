from utils import config, parse_midas_data, sample_utils, diversity_utils, core_gene_utils
import os.path, sys, gzip
import numpy

min_coverage = config.min_median_coverage
min_sample_size = 5
allowed_variant_types = set(['1D','2D','3D','4D'])

# ===========================================================================
# Standard header to read in argument information
# ===========================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--cohort", help="Cohort to consider for computing prevalence", default="all")
parser.add_argument("species", help="Name of specific species to run code on", default="all")
args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
cohort = args.cohort
species_name=args.species
good_species_list = [species_name]

# ===========================================================================

private_snv_directory = '%s/private_snvs/%s/' % (config.data_directory, cohort)
intermediate_filename_template = '%s/%s.txt.gz'
os.system('mkdir -p %s' % private_snv_directory)

# Load subject and sample metadata

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sys.stderr.write("Done!\n")

for species_name in good_species_list:

		sys.stderr.write("Loading haploid samples...\n")

		# Only plot samples above a certain depth threshold that are confidently phaseable.
		snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug, prev_cohort=cohort)
		
		if len(snp_samples) < min_sample_size:
				sys.stderr.write("Not enough haploid samples!\n")
				continue
						
		sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))
		
		core_genes_cohort='all'
		
		sys.stderr.write("Loading whitelisted genes...\n")
		core_genes = core_gene_utils.parse_core_genes(species_name, prev_cohort=core_genes_cohort)
		non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name, prev_cohort=core_genes_cohort)
		shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name, prev_cohort=core_genes_cohort)
		sys.stderr.write("Done! %d core genes and %d shared genes and %d non-shared genes\n" % (len(core_genes), len(shared_pangenome_genes), len(non_shared_genes)))
		
		# Analyze SNPs, looping over chunk sizes. 
		# Clunky, but necessary to limit memory usage on cluster

		# Load SNP information for species_name
		sys.stderr.write("Loading SNPs for %s...\n" % species_name)
		
		private_snvs = []
		
		final_line_number = 0
		while final_line_number >= 0:

				sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
				dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes, prev_cohort=cohort)
				snp_samples = dummy_samples
				sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
				
				chunk_private_snvs = diversity_utils.calculate_private_snvs(snp_samples, allele_counts_map, passed_sites_map)
				
				private_snvs.extend(chunk_private_snvs)
		
		if len(private_snvs)>0:
				
				intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)
				
				# Now add records
				output_file = gzip.open(intermediate_filename,"w")
				# Header
				output_file.write("contig, location, gene_name, var_type, host_id\n")
				for contig, location, gene_name, variant_type, host in private_snvs:
				
						record_str_items = [contig, str(location), gene_name, variant_type, host]           
						record_str = ", ".join(record_str_items)
						output_file.write(record_str)
						output_file.write("\n")
				
				output_file.close()

		sys.stderr.write("Done with %s!\n" % species_name) 

sys.stderr.write("Done looping over species!\n")
