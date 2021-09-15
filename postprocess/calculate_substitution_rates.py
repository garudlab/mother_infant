from utils import sample_utils as su, config, parse_midas_data, diversity_utils, gene_diversity_utils, core_gene_utils
import os.path, sys, gzip, os
import numpy

min_coverage = config.min_median_coverage # 20
min_sample_size = 10

# ===========================================================================
# Standard header to read in argument information
# ===========================================================================

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--cohort", help="Cohort to consider for computing prevalence", default="all")
parser.add_argument("species", help="Name of specific species to run code on")
args = parser.parse_args()

# ===========================================================================
# Different possible prevalence cohorts: cohort and directory names
# 
# all			snp_prevalences/all
# infant	snp_prevalences/infant
# hmp			snp_prevalences/hmp
# mother	snp_prevalences/mother
# ===========================================================================

debug = args.debug
cohort = args.cohort
chunk_size = args.chunk_size
species_name = args.species

# ===========================================================================

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_subject_map = su.parse_sample_subject_map()
sys.stderr.write("Done!\n")

# Header for the output file
record_strs = [", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'Num_muts', 'Num_revs', 'Num_mut_opportunities', 'Num_rev_opportunities'])]

sys.stderr.write("Loading haploid samples...\n")

# Only consider samples above certain depth threshold that are confidently phaseable
snp_samples = diversity_utils.calculate_haploid_samples(species_name, prev_cohort=cohort, debug=debug)

if len(snp_samples) < min_sample_size:
	sys.exit("Not enough haploid samples! (There are %i when minimum is %i)\n" % (len(snp_samples), min_sample_size))

sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

# Load core genes, non-shared reference genes and shared pangenome genes
sys.stderr.write("Loading core genes...\n")

core_genes = core_gene_utils.parse_core_genes(species_name)
non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name, prev_cohort=cohort)
shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name, prev_cohort=cohort)

sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
sys.stderr.write("%d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(core genes only...)\n")
pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])

syn_mut_difference_matrix = numpy.array([]) # 4d sites in core genes
syn_mut_opportunity_matrix = numpy.array([])
syn_rev_difference_matrix = numpy.array([]) # 4d sites in core genes
syn_rev_opportunity_matrix = numpy.array([])

non_mut_difference_matrix = numpy.array([]) # 1d sites in core genes
non_mut_opportunity_matrix = numpy.array([])
non_rev_difference_matrix = numpy.array([]) # 1d sites in core genes
non_rev_opportunity_matrix = numpy.array([])

core_mut_difference_matrix = numpy.array([]) # all sites in core genes
core_mut_opportunity_matrix = numpy.array([])
core_rev_difference_matrix = numpy.array([]) # all sites in core genes
core_rev_opportunity_matrix = numpy.array([])

snp_mut_difference_matrix = numpy.array([]) # all sites in all genes
snp_mut_opportunity_matrix = numpy.array([])
snp_rev_difference_matrix = numpy.array([]) # all sites in all genes
snp_rev_opportunity_matrix = numpy.array([])

final_line_number = 0

while final_line_number >= 0:
	
	sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
	dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size, allowed_genes=non_shared_genes, initial_line_number=final_line_number, prev_cohort=cohort)
	sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
	
	# Calculate fixation matrix
	sys.stderr.write("Calculating matrix of snp differences...\n")
	
	# Synonymous (4D)		
	chunk_syn_mut_difference_matrix, chunk_syn_rev_difference_matrix, chunk_syn_mut_opportunity_matrix, chunk_syn_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['4D']))
	
	# Nonsynonymous (1D) 
	chunk_non_mut_difference_matrix, chunk_non_rev_difference_matrix, chunk_non_mut_opportunity_matrix, chunk_non_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['1D']))
	
	# Core (all)
	chunk_core_mut_difference_matrix, chunk_core_rev_difference_matrix, chunk_core_mut_opportunity_matrix, chunk_core_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes)
	
	# All
	chunk_snp_mut_difference_matrix, chunk_snp_rev_difference_matrix, chunk_snp_mut_opportunity_matrix, chunk_snp_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=non_shared_genes)
	
	sys.stderr.write("Done!\n")
	
	if snp_mut_difference_matrix.shape[0]==0:
		snp_mut_difference_matrix = numpy.zeros_like(chunk_snp_mut_difference_matrix)*1.0
		snp_mut_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		snp_rev_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		snp_rev_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		
		syn_mut_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		syn_mut_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		syn_rev_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		syn_rev_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		
		non_mut_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		non_mut_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		non_rev_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		non_rev_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		
		core_mut_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		core_mut_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		core_rev_difference_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
		core_rev_opportunity_matrix = numpy.zeros_like(snp_mut_difference_matrix)*1.0
	
	# Add syn
	syn_mut_difference_matrix += chunk_syn_mut_difference_matrix
	syn_mut_opportunity_matrix += chunk_syn_mut_opportunity_matrix
	syn_rev_difference_matrix += chunk_syn_rev_difference_matrix
	syn_rev_opportunity_matrix += chunk_syn_rev_opportunity_matrix
	
	# Add non
	non_mut_difference_matrix += chunk_non_mut_difference_matrix
	non_mut_opportunity_matrix += chunk_non_mut_opportunity_matrix
	non_rev_difference_matrix += chunk_non_rev_difference_matrix
	non_rev_opportunity_matrix += chunk_non_rev_opportunity_matrix
	
	# Add core
	core_mut_difference_matrix += chunk_core_mut_difference_matrix
	core_mut_opportunity_matrix += chunk_core_mut_opportunity_matrix
	core_rev_difference_matrix += chunk_core_rev_difference_matrix
	core_rev_opportunity_matrix += chunk_core_rev_opportunity_matrix
	
	# Add all
	snp_mut_difference_matrix += chunk_snp_mut_difference_matrix
	snp_mut_opportunity_matrix += chunk_snp_mut_opportunity_matrix
	snp_rev_difference_matrix += chunk_snp_rev_difference_matrix
	snp_rev_opportunity_matrix += chunk_snp_rev_opportunity_matrix
	
	snp_samples = dummy_samples

# Now calculate gene differences
# Load gene coverage information for this species
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages,     gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples, disallowed_genes=shared_pangenome_genes)
sys.stderr.write("Done! Loaded %d genes\n" % len(gene_names))

gene_sample_list = list(gene_samples)
gene_sample_set = set(gene_samples)

# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")

gene_gain_matrix, gene_loss_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_reads_matrix, gene_depth_matrix, marker_coverages)

good_marker_coverages = (marker_coverages>=min_coverage)

gene_gain_matrix = gene_gain_matrix*good_marker_coverages[:,None]*good_marker_coverages[None,:]
gene_loss_matrix = gene_loss_matrix*good_marker_coverages[:,None]*good_marker_coverages[None,:]
gene_opportunity_matrix = gene_opportunity_matrix*good_marker_coverages[:,None]*good_marker_coverages[None,:]

# Add records to output
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = su.calculate_subject_pairs(sample_subject_map, snp_samples)

for idxs in [same_subject_idxs, diff_subject_idxs]:
	for sample_pair_idx in xrange(0,len(idxs[0])):
		# do both order
		for i,j in [(idxs[0][sample_pair_idx], idxs[1][sample_pair_idx]), (idxs[1][sample_pair_idx], idxs[0][sample_pair_idx])]:
			
			sample_i = snp_samples[i]
			sample_j = snp_samples[j]
			
			record_str_items = [species_name, sample_i, sample_j, '4D', str(syn_mut_difference_matrix[i,j]), str(syn_rev_difference_matrix[i,j]),  str(syn_mut_opportunity_matrix[i,j]), str(syn_rev_opportunity_matrix[i,j])]
			record_strs.append( ", ".join(record_str_items) )
			
			record_str_items = [species_name, sample_i, sample_j, '1D', str(non_mut_difference_matrix[i,j]), str(non_rev_difference_matrix[i,j]),  str(non_mut_opportunity_matrix[i,j]), str(non_rev_opportunity_matrix[i,j])]
			record_strs.append( ", ".join(record_str_items) )
			
			record_str_items = [species_name, sample_i, sample_j, 'core', str(core_mut_difference_matrix[i,j]), str(core_rev_difference_matrix[i,j]),  str(core_mut_opportunity_matrix[i,j]), str(core_rev_opportunity_matrix[i,j])]
			record_strs.append( ", ".join(record_str_items) )
			
			record_str_items = [species_name, sample_i, sample_j, 'all', str(snp_mut_difference_matrix[i,j]), str(snp_rev_difference_matrix[i,j]),  str(snp_mut_opportunity_matrix[i,j]), str(snp_rev_opportunity_matrix[i,j])]
			record_strs.append( ", ".join(record_str_items) )
			
			if (sample_i in gene_sample_set) and (sample_j in gene_sample_set):
					gene_i = gene_sample_list.index(sample_i)
					gene_j = gene_sample_list.index(sample_j)
					
					record_str_items = [species_name, sample_i, sample_j, 'genes', str(gene_loss_matrix[gene_i, gene_j]), str(gene_gain_matrix[gene_i, gene_j]), str(gene_opportunity_matrix[gene_i, gene_j]), '0']
			else:
					# Include samples w/ sufficient SNP coverage that don't have
					# sufficient gene coverage, but mark them with zero opportunities
					record_str_items = [species_name, sample_i, sample_j, 'genes', '0', '0','0','0']
			
			record_strs.append( ", ".join(record_str_items) )

sys.stderr.write("Done with substitution rates for %s!\n" % species_name)

# Write to file

subs_rate_dir = '%s/substitution_rates/%s/' % (config.data_directory, cohort)
os.system('mkdir -p %s' % subs_rate_dir)

sys.stderr.write("Writing intermediate file...\n")
subs_rate_file = '%s/%s.txt.gz' % (subs_rate_dir, species_name)
file = gzip.open(subs_rate_file,"w")
record_str = "\n".join(record_strs)
file.write(record_str)
file.close()
sys.stderr.write("Done!\n")
