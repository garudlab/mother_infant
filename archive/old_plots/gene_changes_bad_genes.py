# Within-host gene change annotation: find bad genes

import sample_utils as su
import config
import parse_midas_data
import os
import pylab
import sys
import numpy
import random

import diversity_utils
import gene_diversity_utils
import sfs_utils
import calculate_substitution_rates
import calculate_temporal_changes
import parse_patric
import species_phylogeny_utils
import core_gene_utils

import stats_utils
from math import log10,ceil
from numpy.random import randint, choice
import pickle

bad_genes = [] # In temporal changes but not in kegg_ids

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species: # empty string evaluates to False
		species_name = other_species
		other_species_str = "_%s" % species_name
else:
		other_species_str = ""
		
################################################################################

if other_species_str=="":
		outFile=open('%sgene_changes_shuffle_all_species.txt' % config.analysis_directory, 'w')
else:
		outFile=open('%sgene_changes_shuffle_%s.txt' % (config.analysis_directory, species_name) , 'w')

modification_difference_threshold = config.modification_difference_threshold # May change
min_coverage = config.min_median_coverage
clade_divergence_threshold = 1e-02 # TODO: change to top level clade definition later
num_trials=100
min_sample_size = 5

within_host_classes = ['gains','losses','all','snps']

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sys.stderr.write("Done!\n")

all_samples = sample_country_map.keys()

# HMP samples
hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

# All mother-infant samples
all_mi_samples = su.get_sample_names('all','all')
mother_samples = su.get_sample_names('mother','all')
infant_samples = su.get_sample_names('infant','all')

# sample set to work with
cohort = infant_samples

if other_species_str == "":
		good_species_list = parse_midas_data.parse_good_species_list()
		if debug:
			good_species_list = [good_species_list[1]]
else:
		good_species_list=[species_name]

# store all the species' data in a dictionary:
all_data={}
#key=species
#value={}, key=gene, valuee=num times gene shows up

for species_name in good_species_list: 
		
		if species_name == 'Escherichia_coli_58110': # Skip to save time, run separately
			continue
		
		dummy_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D']))
		# data structures for storing information for pickling later on
		all_species_gene_changes={}
		#all_species_gene_changes_category={}
		all_species_null={}
		all_data[species_name]={}
		
		####################
		# Analyze the data #
		####################
		
		# Only plot samples above a certain depth threshold that are "haploids"
		haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
		
		# Isolate those samples also belonging to specified set
		haploid_samples = [s for s in haploid_samples if s in cohort]
		
		if len(haploid_samples) < min_sample_size:
				continue
		
		same_subject_idxs = su.calculate_mi_ordered_subject_pairs(sample_order_map, haploid_samples, within_host_type='nonconsecutive', one_per_mi_pair=False)
		# OLD # same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, haploid_samples)
		
		snp_samples = set()
		sample_size = len(same_subject_idxs[0])
		for i, j in zip(same_subject_idxs[0], same_subject_idxs[1]):
				snp_samples.add(haploid_samples[i])
				snp_samples.add(haploid_samples[j])
		
		snp_samples = list(snp_samples)
		allowed_sample_set = set(snp_samples)
		
		if sample_size < min_sample_size:
				continue
		
		# load pre-computed data:
		sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
		substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
		sys.stderr.write("Calculating matrix...\n")
		dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
		snp_samples = dummy_samples
		sys.stderr.write("Done!\n")
		
		sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
		temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
		sys.stderr.write("Done!\n")
		
		# get all genome ids for this species' pan genome:
		genome_ids=parse_midas_data.get_ref_genome_ids(species_name)
		
		# Load the non-shared genes (whitelisted genes):
		non_shared_genes = core_gene_utils.parse_non_shared_pangenome_genes(species_name)
		
		# load the kegg ids for all genomes corresponding to this species:
		kegg_ids=parse_patric.load_kegg_annotations(genome_ids)
		
		##################
		# pangenome null #
		##################
		
		# load all pangenome genes for the species after clustering at 95% identity
		pangenome_gene_names, pangenome_new_species_names=parse_midas_data.load_pangenome_genes(species_name)
		#exclude any genes that are in the whitelisted set from pangenome_gene_names (the pangenome_new_species_names is not used):
		pangenome_gene_names = [gene for gene in pangenome_gene_names if gene not in non_shared_genes]
		
		###########################################
		# load data for between host changes null #
		###########################################
		
		# Load gene coverage information for species_name
		sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
		gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
		sys.stderr.write("Done!\n")
		
		# compute gene cnv for constructing a null based on which genes are present later on.
		gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
		
		# convert gene_samples to list:
		gene_samples=gene_samples.tolist()
		
		# convert gene names to numpy array:
		gene_names=numpy.array(gene_names)
		
		# indexes for different subject pairs
		desired_samples = gene_samples
		
		# Again, only look at samples from specified set
		desired_samples = [s for s in desired_samples if s in cohort]
		
		desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = su.calculate_ordered_subject_pairs(sample_order_map, desired_samples)
		# desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)
		#
		snp_sample_idx_map = su.calculate_sample_idx_map(desired_samples, snp_samples)
		gene_sample_idx_map = su.calculate_sample_idx_map(desired_samples, gene_samples)
		#
		same_subject_snp_idxs = su.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)
		
		#######################
		# within host changes #
		# present gene null -- construct a null consisting of any gene present at either time pt
		#######################
		#
		# store the actual data in this:
		within_host_changes_gene_ids={type:[] for type in within_host_classes}
		#
		# BG: Can't do it this way! Will pick up lots of diploids!
		#for sample_pair in temporal_change_map.keys():
		#		 sample_1=sample_pair[0]
		#		 sample_2=sample_pair[1]
		# 
		for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
				#		 
				i = same_subject_snp_idxs[0][sample_pair_idx]
				j = same_subject_snp_idxs[1][sample_pair_idx]
				#
				sample_i = snp_samples[i] 
				sample_j = snp_samples[j]
				#
				if not ((sample_i in allowed_sample_set) and (sample_j in allowed_sample_set)):
						continue
				#
				# Load SNP changes
				snp_opportunities, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
				#
				# Look at higher threshold if error rate is too high
				if perr>=0.5:
						#
						# Calculate a more fine grained value!
						#
						dfs = numpy.array([0.6,0.7,0.8,0.9])
						perrs = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j,dfs=dfs) * snp_opportunity_matrix[i, j]
						#
						if (perrs<0.5).any():
								# take most permissive one!
								perr_idx = numpy.nonzero(perrs<0.5)[0][0]
								df = dfs[perr_idx]
								perr = perrs[perr_idx]
								#
								# recalculate stuff!		
								perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j,lower_threshold=(1-df)/2.0, upper_threshold=(1+df)/2.0)
						else:
								df = 2
								perr = 1
								mutations = None
								reversions = None
				
				# do same thing for SNP changes
				all_snp_changes = mutations+reversions
				snp_genes = set() # don't double count genes w/ 2 snps. probably same transfer event
				for snp_change in all_snp_changes:
						snp_genes.add(snp_change[0])
				
				# Record the genes missing annotations
				
				for gene in snp_genes:
					if gene not in kegg_ids.keys():
							bad_genes["snp_genes"][species_name].append(gene)
				
				for gene in gene_names:
					if gene not in kegg_ids.keys():
							bad_genes["genes"][species_name].append(gene)				
				
				# Expect pangenome genes to have more since it is 'irrespective of prevalence'?
				for gene in pangenome_gene_names:
					if gene not in kegg_ids.keys():
							bad_genes["pangenome_genes"][species_name].append(gene)

pickle.dump(bad_genes, open("%s/pickles/%s_bad_genes.pkl" % (config.data_directory, species_name), "wb"))
