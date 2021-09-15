from utils import sample_utils, config, parse_midas_data, sfs_utils, diversity_utils, gene_diversity_utils, core_gene_utils
import os, os.path, sys, gzip
import numpy

min_coverage = 20 # config.min_median_coverage
min_sample_size = 2

# ===========================================================================
# Standard header to read in argument information
# ===========================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--min_coverage", type=int, help="min coverage", default=20)
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk_size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--prev_cohort", help="cohort used to compute prevalence", default="all")
parser.add_argument("species", help="Name of specific species to run code on")

args = parser.parse_args()

min_coverage = args.min_coverage
debug = args.debug
chunk_size = args.chunk_size
prev_cohort = args.prev_cohort
species_name=args.species
good_species_list = [species_name]

# ===========================================================================

temporal_change_directory = '%s/temporal_changes_new/cov%i_prev_%s/' % (config.data_directory, min_coverage, prev_cohort)
intermediate_filename_template = '%s/%s.txt.gz'	
os.system('mkdir -p %s' % temporal_change_directory)

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_order_map = sample_utils.parse_sample_order_map()
all_mi_samples = sample_utils.get_sample_names('all','all')

sample_subject_map = sample_utils.parse_sample_subject_map()
mother_samples = sample_utils.get_sample_names('mother','all')
infant_samples = sample_utils.get_sample_names('infant','all')

sys.stderr.write("Done!\n")

intermediate_filename = intermediate_filename_template % (temporal_change_directory, species_name)

output_file = gzip.open(intermediate_filename,"w")
# header!
output_file.write(", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'L','Perr', 'Change1', '...']))
output_file.write("\n")

for species_name in good_species_list:
		
		sample_coverage_map = parse_midas_data.parse_sample_coverage_map(species_name)
		
		sys.stderr.write("Loading SFSs for %s...\t" % species_name)
		samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D']), prev_cohort=prev_cohort)
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading temporal samples...\n")
		# Only plot samples above a certain depth threshold that are involved in timecourse
		snp_samples = diversity_utils.calculate_temporal_samples(species_name)

		# On purpose looking at non-consecutive pairs too
		# (restriction to consecutive pairs is later)
		same_subject_idxs = sample_utils.calculate_mi_ordered_same_subject_pairs(sample_order_map, snp_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)
		
		if len(same_subject_idxs[0]) < min_sample_size:
			sys.exit("Not enough temporal samples! (There are %i when minimum is %i)\n" % (len(same_subject_idxs[0]), min_sample_size))				
		
		sys.stderr.write("Proceeding with %d comparisons of %d temporal samples!\n" % (len(same_subject_idxs[0]), len(snp_samples)))

		# Analyze SNPs, looping over chunk sizes. 
		# Clunky, but necessary to limit memory usage on cluster
		# note core gene prevalence cohort defaults to 'all'

		sys.stderr.write("Loading whitelisted genes...\n")
		non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
		shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
		sys.stderr.write("Done! %d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))

		# Now calculate gene differences
		# Load gene coverage information for species_name
		sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
		gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples, disallowed_genes=shared_pangenome_genes)
		sys.stderr.write("Done!\n")
		
		# Is this redundant?
		same_subject_idxs = sample_utils.calculate_mi_ordered_same_subject_pairs(sample_order_map, snp_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)

		if len(same_subject_idxs[0]) < min_sample_size:
				sys.stderr.write("Not enough temporal samples!\n")
				continue
		
		# Use if HMP # TODO: huh?
		from utils import snps_utils
		private_snv_map = snps_utils.load_private_snv_map(species_name, prev_cohort=prev_cohort)

		# If other dataset, use this (and uncomment private snv line below)
		#import calculate_snp_prevalences
		#snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)				
		
		# Load SNP information for species_name
		sys.stderr.write("Loading SNPs for %s...\n" % species_name)		 
		snp_changes = {}
		gene_changes = {}
		tracked_private_snps = {}
		snp_opportunities = {}
		gene_opportunities = {}
		tracked_private_snp_opportunities = {}
		
		snp_perrs = {}
		gene_perrs = {}
		tracked_private_snp_perrs = {}
		
		snp_difference_matrix = numpy.array([]) # all sites in all genes
		snp_opportunity_matrix = numpy.array([])

		final_line_number = 0

		while final_line_number >= 0:

				sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
				dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number,allowed_genes=non_shared_genes, prev_cohort=prev_cohort)
				sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
				
				# Calculate private snvs
				# (use this if applying to a different dataset)
				# private_snv_map = calculate_private_snvs.calculate_private_snv_map_from_external_snv_prevalence_map(allele_counts_map, snv_freq_map)
		
				# All
				chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map)	 # 

				if snp_difference_matrix.shape[0]==0:
						snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
						snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
						
				# Add all
				snp_difference_matrix += chunk_snp_difference_matrix
				snp_opportunity_matrix += chunk_snp_opportunity_matrix
				
				#same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, snp_samples)
				
				for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):

						i = same_subject_idxs[0][sample_pair_idx]
						j = same_subject_idxs[1][sample_pair_idx]

						sample_i = snp_samples[i]
						sample_j = snp_samples[j]
						
						if not sample_utils.is_same_mi_subject(sample_i, sample_j, sample_subject_map):
							continue
						
						avg_depth_i = sample_coverage_map[sample_i]
						avg_depth_j = sample_coverage_map[sample_j]
						
						chunk_tracked_private_snps = diversity_utils.calculate_tracked_private_snvs(i, j, allele_counts_map, passed_sites_map, avg_depth_i, avg_depth_j, private_snv_map)
						
						chunk_snp_changes = diversity_utils.calculate_snp_differences_between(i, j, allele_counts_map, passed_sites_map, avg_depth_i, avg_depth_j, min_coverage)
		
						sample_pair = (sample_i, sample_j)
		
						if sample_pair not in snp_changes:
								snp_changes[sample_pair] = []
								gene_changes[sample_pair] = []
								snp_opportunities[sample_pair] = 0
								gene_opportunities[sample_pair] = 0
								snp_perrs[sample_pair] = -1
								gene_perrs[sample_pair] = -1
								
								tracked_private_snps[sample_pair] = []
								tracked_private_snp_opportunities[sample_pair] = 0
								tracked_private_snp_perrs[sample_pair] = -1
								
						snp_changes[sample_pair].extend(chunk_snp_changes)
						snp_opportunities[sample_pair] += chunk_snp_opportunity_matrix[i,j]
						
						tracked_private_snps[sample_pair].extend( chunk_tracked_private_snps)
						
						tracked_private_snp_opportunities[sample_pair] += len(chunk_tracked_private_snps)
		
		# Calculate SNP error rate
		for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
				
				i = same_subject_idxs[0][sample_pair_idx]
				j = same_subject_idxs[1][sample_pair_idx]
				
				if not (marker_coverages[i]>=min_coverage)*(marker_coverages[j]>=min_coverage):
						continue
				
				sample_i = snp_samples[i]
				sample_j = snp_samples[j]
				sample_pair = (sample_i, sample_j)
		
				perr = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j)[0]
				
				snp_perrs[sample_pair] = perr
				tracked_private_snp_perrs[sample_pair] = perr

				gene_changes[sample_pair].extend( gene_diversity_utils.calculate_gene_differences_between(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages) )
				
				gene_perr = gene_diversity_utils.calculate_gene_error_rate(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages)[0]
				
				gene_opportunities[sample_pair] = gene_depth_matrix.shape[0]
				
				gene_perrs[sample_pair] = gene_perr

		sys.stderr.write("Done!\n") 

		for sample_i, sample_j in snp_changes.keys():
		
				# First output SNPs
				snp_strs = []
				for snp_change in snp_changes[(sample_i, sample_j)]:				
						
						gene_name, location, variant_type, allele_counts_1, allele_counts_2 = snp_change
						contig = location[0]
						position = location[1]
				
						A1,D1 = allele_counts_1
						A2,D2 = allele_counts_2
				
						snp_str = ('%s;%s;%d;%s;%d;%d;%d;%d' % (gene_name, contig, position, variant_type, A1, D1, A2, D2))
				
						snp_strs.append(snp_str)
				
				record_str_items = [species_name, sample_i, sample_j, 'snps', "%g" % snp_opportunities[(sample_i, sample_j)], "%g" % snp_perrs[(sample_i, sample_j)]] + snp_strs
				record_str = ", ".join(record_str_items)
				output_file.write(record_str)
				output_file.write("\n")
				
				# Now output genes
				gene_strs = []
				for gene_change in gene_changes[(sample_i, sample_j)]:
						gene_idx, coverages_1, coverages_2 = gene_change
						gene_name = gene_names[gene_idx]
						D1,Dm1 = coverages_1
						D2,Dm2 = coverages_2
				
						gene_str = ('%s;%0.2f;%0.2f;%0.2f;%0.2f' % (gene_name, D1, Dm1, D2, Dm2))
						gene_strs.append(gene_str)
		
				record_str_items = [species_name, sample_i, sample_j, 'genes', "%g" % gene_opportunities[(sample_i, sample_j)], "%g" % gene_perrs[(sample_i, sample_j)]] + gene_strs
				record_str = ", ".join(record_str_items)
				output_file.write(record_str)
				output_file.write("\n")

				# Now output private SNPS
				private_snp_strs = []
				for snp_change in tracked_private_snps[(sample_i, sample_j)]:
		
				
						gene_name, location, variant_type, allele_counts_1, allele_counts_2 = snp_change
						contig = location[0]
						position = location[1]
				
						A1,D1 = allele_counts_1
						A2,D2 = allele_counts_2
				
						snp_str = ('%s;%s;%d;%s;%d;%d;%d;%d' % (gene_name, contig, position, variant_type, A1, D1, A2, D2))
				
						private_snp_strs.append(snp_str)
				
				record_str_items = [species_name, sample_i, sample_j, 'private_snps', "%g" % tracked_private_snp_opportunities[(sample_i, sample_j)], "%g" % tracked_private_snp_perrs[(sample_i, sample_j)]] + private_snp_strs
				record_str = ", ".join(record_str_items)
				output_file.write(record_str)
				output_file.write("\n")
		
		sys.stderr.write("Done with %s!\n" % species_name) 

sys.stderr.write("Done looping over species!\n")
output_file.close()
sys.stderr.write("Done!\n")


