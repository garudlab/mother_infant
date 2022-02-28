# Within-host gene change annotation

import sample_utils as su, config, parse_midas_data
import pylab, sys, numpy, random
import diversity_utils, gene_diversity_utils, sfs_utils, calculate_substitution_rates, calculate_temporal_changes, parse_patric, species_phylogeny_utils, core_gene_utils, stats_utils
from numpy.random import choice
import pickle

bad_genes = [] # In temporal changes but not in kegg_ids

#######################################################
# Standard header to read in argument information
#######################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--species', type=str, help='Run the script for one specified species')
parser.add_argument('--sweep-type', type=str, help="Full or partial sweep", default="partial")

args = parser.parse_args()
debug = args.debug
chunk_size = args.chunk_size
species = args.species

sweep_type = args.sweep_type
if sweep_type not in ['full', 'partial']:
	sys.exit("Invalid sweep-type. Choose from full, partial")
if sweep_type == 'full':
	lower_threshold, upper_threshold = 0.2, 0.8
elif sweep_type == 'partial':
	lower_threshold, upper_threshold = 0.35, 0.65
#######################################################

modification_difference_threshold = config.modification_difference_threshold # May change
min_coverage = config.min_median_coverage
clade_divergence_threshold = 1e-02 # TODO: change to top level clade definition later
num_trials=100
min_sample_size = 5

within_host_classes = ['gains','losses','all','snps']

# Gene names for further investigation
genes_of_interest = {gene_name: {} for gene_name in ['Putative large secreted protein SCO0341', 'Uncharacterized sugar:proton symporter', 'Transcriptional regulator, AraC family', 'Transcriptional regulator', 'Two-component system sensor histidine kinase', 'Argininosuccinate lyase (EC 4.3.2.1)', 'Transcriptional regulator YeiE, LysR family', 'Transposase', 'Na+-driven multidrug efflux pump', 'Uncharacterized MFS-type transporter', 'putative membrane protein', 'Type I restriction-modification system, specificity subunit S', 'Acidobacterial duplicated orphan permease (function unknown)', 'putative lipoprotein', 'RNA polymerase ECF-type sigma factor', 'Excinuclease ABC subunit A', 'beta-glycosyl hydrolase', 'Oxidoreductase, aldo/keto reductase family', 'Dienelactone hydrolase and related enzymes', 'Cell surface glycan-binding lipoprotein, utilization system for glycans and polysaccharides (PUL), SusD family', 'Alpha-1,2-mannosidase', 'Integrase', 'Multi antimicrobial extrusion protein (Na(+)/drug antiporter), MATE family of MDR efflux pumps', 'hypothetical protein', 'Uncharacterized protein YphG, TPR-domain containing', 'Outer membrane TonB-dependent transporter, utilization system for glycans and polysaccharides (PUL), SusC family', 'Peptidase M28', 'UDP-glucose 6-dehydrogenase (EC 1.1.1.22)', 'Periplasmic ligand-binding sensor domain COG3292 / BaeS-type histidine kinase / OmpR-type DNA-binding response regulator', 'MSM (multiple sugar metabolism) operon regulatory protein', 'beta-galactosidase (EC 3.2.1.23)', 'Cobalt-zinc-cadmium resistance protein', 'Cystathionine beta-lyase (EC 4.4.1.8)', 'Peptidyl-prolyl cis-trans isomerase (EC 5.2.1.8)']}

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
cohort = hmp_samples + infant_samples

if species is None:
		good_species_list = parse_midas_data.parse_good_species_list()
		if debug:
			good_species_list = [good_species_list[1]]
else:
		good_species_list=[species]

# store all the species' data in a dictionary:
all_data={}
#key=species
#value={}, key=gene, valuee=num times gene shows up

for species_name in good_species_list: 
		
		if species_name == 'Escherichia_coli_58110': # MemoryError... TODO
			continue
		
		try:
			dummy_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D']))
		except: # TODO (Paraprevotella_clara_33712)
			continue
		
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
		try:
			sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
			substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
			sys.stderr.write("Calculating matrix...\n")
			dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
			snp_samples = dummy_samples
			sys.stderr.write("Done!\n")
		except: # TODO Streptococcus_vestibularis_56030
			continue
		
		sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
		temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
		sys.stderr.write("Done!\n")
		
		snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
		sys.stderr.write("Done!\n")
		
		# get all genome ids for this species' pan genome:
		genome_ids=parse_midas_data.get_ref_genome_ids(species_name)
		
		# Load the non-shared genes (whitelisted genes):
		non_shared_genes = core_gene_utils.parse_non_shared_pangenome_genes(species_name)
		
		# load the gene descriptions for all genomes coresponding to this species:
		gene_descriptions=parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
		
		# create gene categories (poor proxy for GO terms):
		gene_categories, gene_category_map = parse_patric.cluster_patric_gene_descriptions(gene_descriptions)
		
		# load the kegg ids for all genomes corresponding to this species:
		kegg_ids=parse_patric.load_kegg_annotations(genome_ids)
		
		# store null data in this to see how the actual data compares. 
		between_host_changes_gene_ids_null={} #dictionary which stores different trials (trial=key)
		present_gene_null={}
		pangenome_null={}
		for change_type in within_host_classes:
				between_host_changes_gene_ids_null[change_type]={}
				present_gene_null[change_type]={}
				pangenome_null[change_type]={}
				for trial in range(0,num_trials):
						between_host_changes_gene_ids_null[change_type][trial]=[]
						present_gene_null[change_type][trial]=[]				
						pangenome_null[change_type][trial]=[]
		
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
		same_subject_gene_idxs = su.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)	 
		#
		diff_subject_snp_idxs = su.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)	 
		diff_subject_gene_idxs = su.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)	 
		#
		between_host_gene_idxs = [] # store idxs of genes that change between hosts
		for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
				snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
				snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
				#
				i = diff_subject_gene_idxs[0][sample_pair_idx]
				j = diff_subject_gene_idxs[1][sample_pair_idx]
				if (marker_coverages[i]>min_coverage) and (marker_coverages[j]>min_coverage):
						if snp_substitution_rate[snp_i, snp_j] < clade_divergence_threshold:
								gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(i,j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
								between_host_gene_idxs.extend(gene_idxs) # collect all gene changes occurring between hosts. Use this for the null.
		#
		#
		#
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
				# Load SNP and gene changes!
				#
				# First SNP changes
				snp_opportunities, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j, lower_threshold, upper_threshold) # Return information on partial vs complete TODO
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
						#		 
						else:
								df = 2
								perr = 1
								mutations = None
								reversions = None
								#
				if mutations==None or perr>=0.5:
						num_mutations = 0
						num_reversions = 0
						num_snp_changes = -1
				else:
						num_mutations = len(mutations)
						num_reversions = len(reversions)
						num_snp_changes = num_mutations+num_reversions
						#
				# Now do gene changes
				gene_opportunities, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
				all_changes=gains+losses
				#
				if (gains==None) or (gene_perr<-0.5) or (gene_perr>0.5):
						num_gains = 0
						num_losses = 0
						num_gene_changes = -1
				else:
						num_gains = len(gains)
						num_losses = len(losses)
						num_gene_changes = num_gains+num_losses
						#
				# Don't want to look at non-modifications or things with high error rates!
				if num_snp_changes<0 or num_snp_changes>=modification_difference_threshold:
						continue
				if (num_snp_changes<=0) and (num_gene_changes<=0):
						continue
				if num_snp_changes <0.5 and num_gene_changes <0.5:
						continue
				gene_change_dictionary={'gains':gains, 'losses':losses, 'all':all_changes}
				#								 
				#iterate through all_changes to store the gene_ids.
				for change_type in ['gains','losses','all']:
						for i in range(0, len(gene_change_dictionary[change_type])):
								gene_id = gene_change_dictionary[change_type][i][0]
								within_host_changes_gene_ids[change_type].append(gene_id)
								
								# Save information on genes of interest
								if gene_id in gene_descriptions:
									gene_name = gene_descriptions[gene_id]
									if gene_name in genes_of_interest:
										sample_pair = frozenset((sample_i, sample_j))
										if sample_pair not in genes_of_interest[gene_name]:
											genes_of_interest[gene_name][sample_pair] = []
										genes_of_interest[gene_name][sample_pair].append(gene_change_dictionary[change_type][i])
				
				# do same thing for SNP changes
				all_snp_changes = mutations+reversions
				for snp_change in all_snp_changes:
						gene_id = snp_change[0]
						# Save information on genes of interest
						if gene_id in gene_descriptions:
							gene_name = gene_descriptions[gene_id]
							if gene_name in genes_of_interest:
								sample_pair = frozenset((sample_i, sample_j))
								if sample_pair not in genes_of_interest[gene_name]:
									genes_of_interest[gene_name][sample_pair] = []
								genes_of_interest[gene_name][sample_pair].append(snp_change)
	
		pickle.dump(genes_of_interest, open("%s/pickles/gene_changes_of_interest_partial.pkl" % config.data_directory, "wb"))
