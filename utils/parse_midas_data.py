import config, sample_utils as su, stats_utils, gene_diversity_utils, midas_db_utils
import sys, bz2, gzip, os.path
import numpy
from math import floor, ceil
from collections import defaultdict

# ==========================================================================
# Set up default source and output directories, shorthand
# ==========================================================================

data_dir = config.data_directory
midas_dir = config.midas_directory
metadata_dir = config.metadata_directory

# We use this one to debug because it was the first one we looked at
debug_species_name = config.debug_species_name

# ==========================================================================
# METHODS FOR PARSING SPECIES METADATA
#		- parse_species_list
#		- parse_depth_sorted_species_list
#		- parse_good_species_list
#		- load_pickled_good_species_list
# ==========================================================================

# ==========================================================================
# Returns a list of all species that MIDAS called SNPS for
# ==========================================================================

def parse_species_list():
	with open("%s/snps/species_snps.txt" % data_dir, 'r') as file:
		species_names = [line.strip() for line in file]
	return species_names

# ==========================================================================
# Returns a list of all species that MIDAS called SNPS for, sorted in order
# of decreasing total sequencing depth [based on marker genes]
# ==========================================================================

def parse_depth_sorted_species_list():
	return parse_global_marker_gene_coverages()[2]

# ==========================================================================
# Returns a list of all species that MIDAS called SNPS for
# that passed a certain depth / prevalence requirement,
# again sorted in order of decreasing total sequencing depth
# ==========================================================================

def parse_good_species_list(min_marker_coverage = config.good_species_min_coverage,
														min_prevalence = config.good_species_min_prevalence):
	good_species_list = []
	
	species_coverage_matrix, samples, species = parse_global_marker_gene_coverages()
	
	for i in range(len(species)):
		species_coverages = species_coverage_matrix[i,:]
		# Number of samples whose marker gene coverage exceeds a threshold
		# must exceed a prevalence threshold
		# Here: at least 10 samples have marker gene coverage >= 10
		if (species_coverages>=min_marker_coverage).sum() >= min_prevalence:
			good_species_list.append(species[i])
	
	return good_species_list

# ==========================================================================
# Returns good_species_list as defined above, but from pickle file
# ==========================================================================

def load_pickled_good_species_list():
	
	import pickle	
	pickle_path = '%s/pickles/good_species_list.pkl' % data_dir
	
	if os.path.isfile(pickle_path):
		return pickle.load(open(pickle_path, 'rb'))
	else:
		good_species_list = parse_good_species_list()
		pickle.dump(good_species_list, open(pickle_path, 'wb'))
		return good_species_list

# ==========================================================================
# METHODS FOR PARSING COVERAGE OF DIFFERENT SPECIES ACROSS SAMPLES
#		- parse_global_marker_gene_coverages()
#		- parse_species_marker_gene_coverages(species)
#		- parse_marker_gene_coverage_distribution(species)
#		- parse_gene_coverages(species)
#		- parse_coverage_distribution(species)
#		- parse_sample_coverage_map(species)
# ==========================================================================

# ==========================================================================
# Loads marker gene coverages produced by MIDAS
# for all species in which SNPs were called 
#
# Returns: species-by-sample matrix of marker gene coverages,
#					 with species sorted in descending order of total coverage;
#					 ordered list of sample ids; ordered list of species names
# ==========================================================================

def parse_global_marker_gene_coverages():
	
	desired_species_names = set(parse_species_list())
	
	# coverage.txt: average read-depth of 15 marker genes per species
	# (total bp of mapped reads/total bp of 15 marker-genes)
	
	file = bz2.BZ2File("%s/species/coverage.txt.bz2" % data_dir, 'r')
	samples = file.readline().strip().split()[1:] # header
	samples = su.parse_merged_sample_names(samples) # remove c's
	
	species = [] # list of species
	species_coverage_matrix = [] # rows - species, cols - samples
	
	for line in file:
		items = line.strip().split()
		species_name = items[0]
		coverages = numpy.array([float(cov) for cov in items[1:]])
		
		if species_name in desired_species_names:
			species.append(species_name)
			species_coverage_matrix.append(coverages)
	
	file.close()
	
	# Sort by marker gene coverage, summed across all samples, from high to low	
	species, species_coverage_matrix = zip(*sorted(zip(species, species_coverage_matrix), key=lambda pair: pair[1].sum(), reverse=True))
	
	species_coverage_matrix = numpy.array(species_coverage_matrix)
	
	return species_coverage_matrix, samples, species

# ==========================================================================
# Loads marker gene coverages produced by MIDAS for a particular species
# Returns: list of average marker coverages, list of samples
# ==========================================================================

def parse_species_marker_gene_coverages(desired_species_name):
	
	species_coverage_matrix, samples, species = parse_global_marker_gene_coverages()
	
	try:
		species_idx = species.index(desired_species_name)
		return species_coverage_matrix[species_idx,:], samples
	except:
		return None # desired species name not found

# ==========================================================================
# Loads site coverages within marker genes for a particular species
# Returns dictionary: sample -> marker gene -> (location, coverage)
# 
# Note that marker_coverage_distribution.txt
# is from calculate_marker_coverage_distribution
# ==========================================================================

def parse_marker_gene_coverage_distribution(desired_species_name):
	
	coverage_distribution_file = bz2.BZ2File("%s/snps/%s/marker_coverage_distribution.txt.bz2" % (data_dir, desired_species_name))
	coverage_distribution_file.readline() # header
	
	# Dictionary: sample -> marker gene -> (location, coverage) list
	marker_cov_dict = defaultdict(dict)
	
	for line in coverage_distribution_file:
		items = line.strip().split("\t")
		sample, gene_name = items[0].split(",")
		sample = sample[:-1] if (sample[-1] == 'c') else sample # remove c
		
		locations, coverages = [], []
		
		for item in items[1:]:
			loc, cov = item.split(",")
			locations.append(long(loc))
			coverages.append(float(cov))
		
		marker_cov_dict[sample][gene_name] = (numpy.array(locations), numpy.array(coverages))
	
	return marker_cov_dict

# ==========================================================================
# Loads gene coverages for a particular species (gene_name -> avg coverage)
# 
# Recall how gene_coverage.txt was generated:
# Average depth, i.e. sum of depths across all sites with known variant
# type divided by the number of such sites, is reported for each gene
# in each sample (so long some sample-gene has at least one such site).
# ==========================================================================

def parse_gene_coverages(desired_species_name):

	gene_cov_file = bz2.BZ2File("%s/snps/%s/gene_coverage.txt.bz2" % (data_dir, desired_species_name))

	samples = gene_cov_file.readline().strip().split()[1:] # header
	samples = su.parse_merged_sample_names(samples) # remove c's
	
	# Dictionary: pangenome gene -> average depth in each sample
	gene_coverages = {}
	
	for line in gene_cov_file:
		items = line.split()
		gene_name = items[0]
		depths = numpy.array([float(item) for item in items[1:]])
		gene_coverages[gene_name] = depths
	
	return gene_coverages, samples

# ==========================================================================
# Recall how the site coverage distribution files were generated:
# 
#		full_coverage_distribution.txt: for each sample, consider all sites
#		with known variant type (1D, 2D, 3D, or 4D) and output number of sites
#		which have each depth genome-wide
# 
#		coverage_distribution.txt: consider only sites with known variant type
#		and that have coverage >= 3 in at least 95% of all samples
# 		cohort options: infant, adult, hmp, any
# 
# This function loads coverage histogram (dictionary: depth -> site count)
# for each sample for a given species.
# ==========================================================================

def parse_coverage_distribution(desired_species_name, prevalence_filter = True, remove_c = True, prev_cohort='any'):
	
	if prev_cohort in ['infant', 'adult', 'hmp']:
		prev_cohort_str = 'prev_' + prev_cohort
	else:
		prev_cohort_str = prev_cohort
	
	full_str = "" if prevalence_filter else "full_"
	cov_dist_file = bz2.BZ2File("%s/snps/%s/%scoverage_distribution_%s.txt.bz2" % (data_dir, desired_species_name, full_str, prev_cohort_str))

	cov_dist_file.readline() # remove header
	samples, sample_coverage_histograms = [], []
	for line in cov_dist_file:
		items = line.strip().split()
		sample_coverage_histogram = {}
		for item in items[1:]:
			subitems = item.split(",")
			sample_coverage_histogram[float(subitems[0])] = float(subitems[1])
		sample_coverage_histograms.append(sample_coverage_histogram)
		samples.append(items[0])
	
	if remove_c == True:
		samples = su.parse_merged_sample_names(samples)		 
	
	return sample_coverage_histograms, samples

# ==========================================================================
# Returns median site coverage (i.e. half of all sites have same or lower)
# for each sample in a dictionary (sample -> median_coverage)
# ==========================================================================

def parse_sample_coverage_map(species, prev_cohort='any'):
	
	sample_coverage_histograms, samples = parse_coverage_distribution(species, prev_cohort=prev_cohort)
	
	sample_coverage_map = {}
	
	for hist, sample in zip(sample_coverage_histograms, samples):
		median_coverage = stats_utils.calculate_nonzero_median_from_histogram(hist)
		sample_coverage_map[sample] = median_coverage
	
	return sample_coverage_map

# ==========================================================================
# METHODS FOR LOADING SETS OF GENES FOR A SPECIES
# (technically, these should be in midas_db_utils because they do not
# depend on a particular dataset. But I'll keep these here anyways...)
# 
#		- load_pangenome_genes
#		- load_reference_genes
#		- load_marker_genes
#		- load_core_genes
# ==========================================================================

# ==========================================================================
# Loads MIDAS pangenome (after clustering at 95% identity) for given species
# Returns a complete set of genes (irrespective of prevalence)
# ==========================================================================

def load_pangenome_genes(species_name):
	
	# Open MIDAS output: presence absence calls for all pangenome centroids
	gene_presabs_file =	 bz2.BZ2File("%s/genes/%s/genes_presabs.txt.bz2" % (data_dir, species_name), 'r')
	gene_presabs_file.readline() # remove header	
	genes = [line.split()[0] for line in gene_presabs_file]
	
	# Maps original centroid to reference-corrected centroid
	centroid_gene_map = midas_db_utils.load_centroid_gene_map(species_name)
	ref_corrected_genes = [centroid_gene_map[c] for c in genes]
	
	return set(genes), set(ref_corrected_genes)

# ==========================================================================
# Loads list of genes in reference genome used by MIDAS for a given species
# ==========================================================================

def load_reference_genes(desired_species_name):
	
	features_file = gzip.open("%s/rep_genomes/%s/genome.features.gz" % (midas_dir, desired_species_name), 'r')	
	features_file.readline() # remove header
	reference_genes = [line.split()[0] for line in features_file]
	return set(reference_genes)

# ===========================================================================
# Loads list of MIDAS marker genes for a given species
# (by default, load only those 15 IDs that are in the reference genome)
# ===========================================================================

def load_marker_genes(desired_species_name, require_in_reference_genome=True):	 
	
	# Chosen markers (From table S7 of Nayfach et al (Genome Res 2016)
	marker_ids = set(['B000032', 'B000039', 'B000041', 'B000062', 'B000063', 'B000065', 'B000071', 'B000079', 'B000080', 'B000081', 'B000082', 'B000086', 'B000096', 'B000103', 'B000114'])
	
	# Get list of reference genome genes
	reference_genes = set(load_reference_genes(desired_species_name))
	
	# Maps marker ID to gene name
	marker_gene_map = {marker_id : [] for marker_id in marker_ids}
	marker_genes = []
	
	# MIDAS marker gene database: gene ID, rep(?) genome ID, marker ID
	# Note that multiple gene IDs may correspond to the same marker gene
	marker_gene_file = open("%s/marker_genes/phyeco.map" % (midas_dir), 'r')
	marker_gene_file.readline() # header
	
	for line in marker_gene_file:
		gene_name, gene_len, genome_id, species_name, marker_id = line.strip().split("\t")
		
		# If you only want marker genes in the reference genome,
		# the gene name must be in reference genes
		if species_name == desired_species_name and marker_id in marker_ids and ((not require_in_reference_genome) or (gene_name in reference_genes)):
			marker_genes.append(gene_name)
			marker_gene_map[marker_id].append(gene_name)
	
	return set(marker_genes)
		 
# ===========================================================================
# Loads a subset of "core" genes
# prev_cohort: adult, hmp, infant, all, any
# *Deprecated: Use core_gene_utils.parse_core_genes instead
# ===========================================================================	

def load_core_genes(desired_species_name, prev_cohort='all'):
	import core_gene_utils
	return core_gene_utils.parse_core_genes(desired_species_name, prev_cohort)

# ==========================================================================
# METHODS THAT USE MIDAS INTERMEDIATE FILES 
# ==========================================================================

# ==========================================================================
# Returns read counts of allowed genes for each genome
# in each sample-species pair; uses MIDAS intermediate files 
# ==========================================================================

def parse_99_percent_genes(species, samples, allowed_genes=[]):

	# Dictionary: sample -> gene -> read count
	data = defaultdict(dict)
	
	# Assume samples are unmerged
	for sample in samples:
		int_genes_output_dir = "%s/%s/genes/output" % (config.int_data_directory, sample)
		try:
			file = gzip.open("%s/%s.genes.gz" % (int_genes_output_dir, species), 'r')
			file.readline() # remove header
			for line in file:
				gene, count_reads, coverage, copynum = line.strip().split('\t')
				if gene in allowed_genes:
					data[sample][gene] = int(count_reads)
		except:
			print("%s has no genes output for sample %s" % (species, sample))
	
	# Dictionary: gene -> list of read counts, ordered by sample
	data_numpy_array_dict = defaultdict(list)
	
	for gene in allowed_genes:
		for sample in samples:
			read_count = data[sample][gene] if (gene in data[sample]) else 0
			data_numpy_array_dict[gene].append(read_count)
		data_numpy_array_dict[gene] = numpy.asarray(data_numpy_array_dict[gene])
	
	# Dictionary: ref_genome -> list of read counts, ordered by gene, sample
	ref_genome_dict = defaultdict(list)
	for gene in allowed_genes:
		ref_genome='.'.join(gene.split('.')[0:2])
		ref_genome_dict[ref_genome] += data_numpy_array_dict[gene]
	
	return ref_genome_dict

# ==========================================================================
# Parse MIDAS intermediate species file to get a list of species
# at least 3x coverage
# ==========================================================================

def parse_intermediate_species_file(sample_id, inFN):
	inFile = open(inFN,'r')
	inFile.readline() # remove header
	
	species_list=[]
	
	for line in inFile:
		items=line.strip().split('\t')
		species_id=items[0]
		coverage=float(items[2])
		if coverage >=3.0: 
			species_list.append(species_id)
	
	return set(species_list)

# ==========================================================================
# COMPLICATED POSTPROCESSING METHODS...
# ==========================================================================

def calculate_relative_depth_threshold_map(sample_coverage_histograms, samples, min_nonzero_median_coverage=5, lower_factor=0.5, upper_factor=2):
	
	# returns map of sample name: coverage threshold
	# essentially filtering out samples whose marker depth coverage
	# does not exceed the average coverage threshold
	
	depth_threshold_map = {}
	for i in xrange(0,len(samples)):
		
		# Check if coverage distribution meets certain requirements
		is_bad_coverage_distribution = False
		
		# First check if passes median coverage requirement
		nonzero_median_coverage = stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histograms[i])
		if round(nonzero_median_coverage) < min_nonzero_median_coverage:
				is_bad_coverage_distribution=True
		
		# Passed median coverage requirement
		# Now check whether a significant number of sites fall between lower and upper factor. 
		lower_depth_threshold = floor(nonzero_median_coverage*lower_factor)-0.5 # why is 0.5 being added/subtracted? NRG
		upper_depth_threshold = ceil(nonzero_median_coverage*upper_factor)+0.5
		
		depths, depth_CDF = stats_utils.calculate_CDF_from_histogram(sample_coverage_histograms[i])
		# remove zeros
		if depths[0]<0.5:
				depth_CDF -= depth_CDF[0]
				depth_CDF /= depth_CDF[-1]
		
		fraction_in_good_range = depth_CDF[(depths>lower_depth_threshold)*(depths<upper_depth_threshold)].sum()
		
		if fraction_in_good_range < 0.6: #where does 0.6 come from? NRG
				is_bad_coverage_distribution=True
		
		if is_bad_coverage_distribution:
				lower_depth_threshold = 1000000001
				upper_depth_threshold = 1000000001
		
		depth_threshold_map[samples[i]] = (lower_depth_threshold, upper_depth_threshold)
	
	return depth_threshold_map

# ===========================================================================
# Reads midas output and prints to stdout in a format 
# suitable for further downstream processing
#
# In the process, filters sites that fail to meet the depth requirements
# ===========================================================================

def pipe_snps(species_name, min_nonzero_median_coverage=config.pipe_snps_min_nonzero_median_coverage, lower_factor=config.pipe_snps_lower_depth_factor, upper_factor=config.pipe_snps_upper_depth_factor, min_samples=config.pipe_snps_min_samples, debug=False):
	
	# lower_factor = 0.3 is the default to be consistent with MIDAS gene presence criterion
	# upper factor = 3 is the default for (logarithmic) symmetry 
	# min_samples=4 is the default because then the site is guaranteed to be present in 
	# at least 2 independent people. 
	# NRG: Why is a site guaranteed to be in at least 2 independent people?
	# BG: In our cohort, the maximum number of samples per person is 3. 
	#			If there are 4 samples, then they must be spread across at least 2 people. 
	
	# Load genomic coverage distributions
	sample_coverage_histograms, sample_list = parse_coverage_distribution(species_name, remove_c=False)
	
	# depth threshold map returns the lower and upper depth values that are 0.3*median and 3*median depth in the data. 
	depth_threshold_map = calculate_relative_depth_threshold_map(sample_coverage_histograms, sample_list, min_nonzero_median_coverage, lower_factor, upper_factor)
	
 
	# Open MIDAS output files
	ref_freq_file = bz2.BZ2File("%s/snps/%s/snps_ref_freq.txt.bz2" % (data_dir, species_name),"r")
	depth_file = bz2.BZ2File("%s/snps/%s/snps_depth.txt.bz2" % (data_dir, species_name),"r")
	alt_allele_file = bz2.BZ2File("%s/snps/%s/snps_alt_allele.txt.bz2" % (data_dir, species_name),"r")
	info_file = bz2.BZ2File("%s/snps/%s/snps_info.txt.bz2" % (data_dir, species_name),"r")
	marker_file = bz2.BZ2File("%s/snps/%s/marker_coverage.txt.bz2" % (data_dir, species_name))
	
	# get header lines from each file
	depth_line = depth_file.readline()
	ref_freq_line = ref_freq_file.readline()
	alt_line = alt_allele_file.readline()
	info_line = info_file.readline()
	marker_line = marker_file.readline()
	
	# get list of samples
	depth_items = depth_line.split()
	samples = numpy.array(depth_items[1:])
	
	# BHG (06/24/17) removed this so that all "raw" data have "c"s. 
	# All functions that write something keep them. 
	# All loading functions strip them. 
	#samples= parse_merged_sample_names(samples) # NRG (06/06/17): I added this so that the keys in dictionary are compatible. 

	# samples
	prevalence_threshold = min([min_samples*1.0/len(samples), 0.5])
	
	# create depth threshold vector from depth threshold map
	lower_depth_threshold_vector = []
	upper_depth_threshold_vector = []
	for sample in samples:
			lower_depth_threshold_vector.append(depth_threshold_map[sample][0])
			upper_depth_threshold_vector.append(depth_threshold_map[sample][1])
			
	lower_depth_threshold_vector = numpy.array(lower_depth_threshold_vector)
	upper_depth_threshold_vector = numpy.array(upper_depth_threshold_vector)
	
	# Figure out which samples passed our avg_depth_threshold
	passed_samples = (lower_depth_threshold_vector<1e09) #1e09 comes from the calculate_relative_depth_threshold_map definition above, which is a code for a bad sample. A bad sample has median depth less than 5 or greater than 0.6 fraction of the genome is outside the acceptable range of good depths. 
	total_passed_samples = passed_samples.sum()
	
	# Let's focus on those from now on
	samples = list(samples[passed_samples])
	lower_depth_threshold_vector = lower_depth_threshold_vector[passed_samples]
	upper_depth_threshold_vector = upper_depth_threshold_vector[passed_samples]
	
	#print lower_depth_threshold_vector
	
	# print header
	print_str = "\t".join(["site_id"]+samples)
	print print_str
	
	# Only going to look at 1D, 2D, 3D, and 4D sites
	# (we will restrict to 1D and 4D downstream)
	allowed_variant_types = set(['1D','2D','3D','4D'])
	
	allele_counts_syn = [] # alt and reference allele counts at 4D synonymous sites with snps
	locations_syn = [] # genomic location of 4D synonymous sites with snps
	genes_syn = [] # gene name of 4D synonymous sites with snps
	passed_sites_syn = numpy.zeros(len(samples))*1.0
	
	allele_counts_non = [] # alt and reference allele counts at 1D nonsynonymous sites with snps
	locations_non = [] # genomic location of 1D nonsynonymous sites
	genes_non = [] # gene name of 1D nonsynonymous sites with snps
	passed_sites_non = numpy.zeros_like(passed_sites_syn)
	
	num_sites_processed = 0
	while True:
					
			# load next lines
			depth_line = depth_file.readline()
			ref_freq_line = ref_freq_file.readline()
			alt_line = alt_allele_file.readline()
			info_line = info_file.readline()
			
			# quit if file has ended
			if depth_line=="":
					break
			
			# parse site info
			info_items = info_line.split("\t")
			variant_type = info_items[5]
			
			# make sure it is either a 1D or 4D site
			if not variant_type in allowed_variant_types:
					continue
	
			# continue parsing site info
			gene_name = info_items[6]
			site_id_items = info_items[0].split("|")
			# NRG: added this if condition to deal with extra 'accn' in db swap. 
			if site_id_items[0]=='accn':
					contig = site_id_items[1]
					location = site_id_items[2]
			else:
					contig = site_id_items[0] 
					location = site_id_items[1]
			
	
			# now parse allele count info
			depths = numpy.array([float(item) for item in depth_line.split()[1:]])
			ref_freqs = numpy.array([float(item) for item in ref_freq_line.split()[1:]])
			
			depths = depths[passed_samples]
			ref_freqs = ref_freqs[passed_samples]
			
			#if (ref_freqs==1.0).all():
			#		 sys.stderr.write("Non-polymorphic site!\n")
			
			refs = numpy.round(ref_freqs*depths)	 
			alts = depths-refs
			
			passed_sites = (depths>=lower_depth_threshold_vector)*1.0
			passed_sites *= (depths<=upper_depth_threshold_vector)
			
			#print passed_sites.sum(), total_passed_samples, passed_sites.sum()/total_passed_samples
			
			# make sure the site is prevalent in enough samples to count as "core"
			if (passed_sites).sum()*1.0/total_passed_samples < prevalence_threshold:
					continue
					#passed_sites *= 0
					
			refs = refs*passed_sites
			alts = alts*passed_sites
			depths = depths*passed_sites
			
			total_alts = alts.sum()
			total_refs = depths.sum()
			total_depths = total_alts+total_refs
			
			# BG: 05/18: moving the polarization part to another part of the pipeline
			# so that we can use HMP polarization with other datasets. 
			# at the moment, still saving polarization state.
			
			polarization = "R"		
			new_site_id_str = "|".join([contig, location, gene_name, variant_type, polarization])
			
			
			# print string
			read_strs = ["%g,%g" % (A,A+R) for A,R in zip(alts, refs)]
			print_str = "\t".join([new_site_id_str]+read_strs)
			
			print print_str
			#print total_alts
			
			num_sites_processed+=1
			if num_sites_processed%10000==0:
					#sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))	 
					if debug:
							break
	
	ref_freq_file.close()
	depth_file.close()
	alt_allele_file.close()
	info_file.close()
	
	# returns nothing

# ===========================================================================
# Loads list of SNPs and counts of target sites from annotated SNPs file
# for a particular species, from allowed samples at allowed sites
#
# returns:	desired_samples - samples considered (default all)
#						allele_counts_map - (gene -> vartype -> locations, allele counts)
#						passed_sites_map - (gene -> vartype -> location -> pass matrix)
#						final_line_number - last line in annotated_snps.txt processed
# ===========================================================================

def parse_snps(species_name, debug=False, allowed_samples=[], allowed_genes=[], allowed_variant_types=['1D','2D','3D','4D'], initial_line_number=0, chunk_size=1000000000, alt_passed_sites_map=False, prev_cohort='all'):
		
		from utils import snps_utils
		
		# Load population freqs (for polarization purposes)		 
		population_freqs = snps_utils.parse_population_freqs(prev_cohort, species_name, polarize_by_consensus=False)
	 
		# Open post-processed MIDAS output
		snps_dir = "%s/snps/%s/" % (data_dir, species_name)
		snp_file = bz2.BZ2File("%s/annotated_snps.txt.bz2" % snps_dir, 'r')
		
		# Get lists of desired samples, genes, and variant types
		items = snp_file.readline().strip().split()[1:]
		samples = su.parse_merged_sample_names(items)
		samples_list = list(samples)
		
		if len(allowed_samples) == 0:
			allowed_sample_set = set(samples)
		else:
			allowed_sample_set = (set(allowed_samples) & set(samples))
		
		desired_sample_idxs = numpy.array(sorted([samples_list.index(sample) for sample in allowed_sample_set]))
		desired_samples =	 samples[desired_sample_idxs]
		
		allowed_genes = set(allowed_genes)
		allowed_variant_types = set(allowed_variant_types)
		
		# Map: gene name -> variant type -> (list of locations, matrix of allele counts)
		allele_counts_map = {}
		
		# Map: gene name -> variant type -> location -> sample-sample matrix of whether both samples can be called at that site
		passed_sites_map = {}
		
		num_sites_processed = 0
		line_number = -1
		final_line_number = -1
		previous_gene_name = ""
		gene_name = ""
		
		for line in snp_file:
				
				line_number += 1 # Start at line 0 (really, second line of file)
				previous_gene_name = gene_name
				
				if line_number < initial_line_number:
						continue # Skip until desired initial line number is reached
				
				items = line.split()
				
				# Load information about site
				info_items = items[0].split("|")
				chromosome = info_items[0]
				location = long(info_items[1])
				gene_name = info_items[2]
				variant_type = info_items[3]
				
				if len(info_items) > 5: # for backwards compatability
						polarization = info_items[4]
						pvalue = float(info_items[5])
				else: 
						polarization="R" # not correct, but avoids a crash (TODO: ??)
						pvalue = float(info_items[4])
				
				if num_sites_processed >= chunk_size and gene_name != previous_gene_name:
						# Chunk processed, we are done for now!
						final_line_number = line_number
						break
				
				# Ignore site if variant type or gene is not allowed
				
				if not variant_type in allowed_variant_types:
						continue
						
				if len(allowed_genes)>0 and (not gene_name in allowed_genes):
						continue
				
				# Load alt and depth counts for all desired samples
				alts, depths = [], []
				
				for idx in desired_sample_idxs:		 
					item = items[1+idx]
					subitems = item.split(",")
					alts.append(float(subitems[0]))
					depths.append(float(subitems[1]))
				
				alts = numpy.array(alts)
				depths = numpy.array(depths)
				
				# Obtain population frequency of alt allele
				if (chromosome, location) in population_freqs:
					population_freq = population_freqs[(chromosome, location)]
				else: # alt population prevalence is (probably? TODO) 0
					population_freq = 0
				
				# Polarize SFS according to population freq
				if population_freq > 0.5: # This means alt allele is the major allele
						alts = depths - alts
						polarization = 'A'
				
				# 0/1 array in which sample has 1 if site has nonzero depth, else 0
				passed_sites = (depths>0) * 1.0
				
				# Set up passed_sites_map, allele_counts_map for new genes
				if gene_name not in passed_sites_map:
					if alt_passed_sites_map == True:
						passed_sites_map[gene_name] = {v: {'locations':[], 'sites':[]} for v in allowed_variant_types}
					else:
						passed_sites_map[gene_name] = {v: {'locations':[], 'sites':numpy.zeros((len(desired_samples), len(desired_samples)))} for v in allowed_variant_types}
					allele_counts_map[gene_name] = {v: {'locations':[], 'alleles':[]} for v in allowed_variant_types}
				
				# Store passed sites information
				full_location = (chromosome,location)
				passed_site_matrix = passed_sites[:,None]*passed_sites[None,:]
				
				passed_sites_map[gene_name][variant_type]['locations'].append(full_location)
				if alt_passed_sites_map == True:
					passed_sites_map[gene_name][variant_type]['sites'].append(passed_site_matrix)
				else:
					passed_sites_map[gene_name][variant_type]['sites'] += passed_site_matrix
				
				# Zero out non-passed sites (TODO: pointless?)	 
				alts = alts*passed_sites
				depths = depths*passed_sites
				
				# Calculate whether SNP has passed, i.e. number of alternate alleles
				# make up more than 5% of all alleles in at least one sample
				# and pvalue < 0.05 (where pvalue is from annotated_snps.txt)
				
				alt_threshold = numpy.ceil(depths*config.parse_snps_min_freq) + 0.5
				snp_passed = ((alts > alt_threshold).sum()>0) and (pvalue<0.05)
				
				# Criteria used in Schloissnig et al (Nature, 2013)
				#total_alts = alts.sum()
				#total_depths = depths.sum()
				#pooled_freq = total_alts/((total_depths+(total_depths==0))
				#snp_passed = (freq>0.01) and (total_alts>=4) and ((total_depths-total_alts)>=4)
				
				# Store allele counts only if the site is interesting
				if snp_passed:
					allele_counts = numpy.transpose(numpy.array([alts,depths-alts]))
					allele_counts_map[gene_name][variant_type]['locations'].append(full_location)
					allele_counts_map[gene_name][variant_type]['alleles'].append(allele_counts)
			
					num_sites_processed += 1
					
					if num_sites_processed>0 and num_sites_processed%1000==0:
							sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))		
							if debug:
									break
		
		snp_file.close()
		
		for gene_name in passed_sites_map.keys():
			for variant_type in passed_sites_map[gene_name].keys():
				allele_counts_map[gene_name][variant_type]['alleles'] = numpy.array(allele_counts_map[gene_name][variant_type]['alleles'])
		
		return desired_samples, allele_counts_map, passed_sites_map, final_line_number

# ===========================================================================
# Calculates within-sample synonymous pi directly from annotated_snps.txt.bz2 
# Ugly hack (since it does not encourage code re-use and puts pop-gen logic
# in the basic parsing scripts) but we need it so that we can call parse_snps 
# on subsets of samples later on to improve performance
#
# Returns: samples, vector pi_s (raw counts), vector of opportunities
# ===========================================================================

def parse_within_sample_pi_new(species_name, allowed_genes=set([]), allowed_variant_types=set(['4D']), debug=False):
	
	samples, sfs_map = parse_within_sample_sfs(species_name, allowed_variant_types)
	
	total_pi = []
	total_opportunities = []
	
	for sample in samples: 
			p,n = diversity_utils.calculate_pi_from_sfs_map(sfs_map[sample])
			total_pi.append(p)
			total_opportunities.append(n)
	
	total_pi = numpy.array(total_pi)*1.0
	total_opportunities = numpy.array(total_opportunities)*1.0
	
	return samples, total_pi, total_opportunities
		
# ===========================================================================
# Calculates within-sample sfs directly from annotated_snps.txt.bz2. 
# Ugly hack (since it does not encourage code re-use and puts pop-gen logic
# in the basic parsing scripts) but we need it so that we can call parse_snps 
# on subsets of haploid samples later on to improve performance
# 
# prev cohorts: all, hmp, infant
# where core genes are still based on all but prev cohort is for snp prevalence
#
# returns vector of samples, vector of sfs maps
# ===========================================================================

def parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D']), prev_cohort='all'):
	
	# First read (filtered) genome-wide coverage distribution
	sfs_file = bz2.BZ2File("%s/snps/%s/within_sample_sfs_prev_%s.txt.bz2" % (data_dir, species_name, prev_cohort), 'r')
	sfs_file.readline() # header
	
	sfs_map = defaultdict(dict)
	samples = []
	
	for line in sfs_file:
		items = line.strip().split("\t")
		sample, variant_type, sfs_items = items[0], items[1], items[2:]
		
		if variant_type not in allowed_variant_types:
			continue
		
		sample = sample[:-1] if sample[-1] == 'c' else sample
		if sample not in sfs_map:
			samples.append(sample)
		
		for sfs_item in sfs_items:
			subitems = sfs_item.split(",")
			D = long(subitems[0]) # Total depth
			A = long(subitems[1]) # Alt allele depth
			n = long(subitems[2]) # Number of sites with this D,A
			reverse_n = float(subitems[3]) # Reverse count (?)
			
			if D<0.5: # Ignore 0,0
				continue
			
			if (A,D) not in sfs_map[sample]: # why is this A,D instead of D,A?
				sfs_map[sample][(D,A)] = [0,0.0]
			
			sfs_map[sample][(D,A)][0] += n
			sfs_map[sample][(D,A)][1] += reverse_n
	
	return numpy.array(samples), sfs_map

# ===========================================================================
# Loads MIDAS's pangenome coverage data for a given species
# returns (lots of things, see below)
# ===========================================================================

def parse_pangenome_data(species_name, allowed_samples = None, allowed_genes=[], convert_centroid_names = True, disallowed_genes= []):
		
		gene_reads_path =	"%s/genes/%s/genes_reads.txt.bz2" % (data_dir, species_name)
		
		if not os.path.isfile(gene_reads_path):
				return [], [], [], [], [], []
				
		# Open post-processed MIDAS output
		
		# genes_reads.txt: number of reads mapped to each gene per sample
		# TODO: is a single base overlap sufficient?
		gene_reads_file =	 bz2.BZ2File("%s/genes/%s/genes_reads.txt.bz2" % (data_dir, species_name),"r")
		
		# genes_depth.txt: average read depth of each gene per sample
		# TODO: averaged over all bases?
		gene_depth_file =	 bz2.BZ2File("%s/genes/%s/genes_depth.txt.bz2" % (data_dir, species_name),"r")
		
		# genes_presabs.txt: presence (1)/absence (0) of each gene per sample
		# where presabs is based on copy number threshold
		gene_presabs_file =	 bz2.BZ2File("%s/genes/%s/genes_presabs.txt.bz2" % (data_dir, species_name),"r")
		
		# Get marker gene coverages
		# marker_covs, mcov_samples = parse_species_marker_gene_coverages(species_name)
		# # Convert to dictionary form
		# marker_coverage_map = {}
		# for mcov, sample in zip(marker_covs, mcov_samples):
		#		marker_coverage_map[sample] = mcov
		# Old way of getting marker gene coverages from genes_summary
		gene_summary_file = file("%s/genes/%s/genes_summary.txt" % (data_dir, species_name),"r")
		marker_coverage_map = {}
		gene_summary_file.readline() # header
		marker_coverage_samples = []
		marker_coverages = []
		for summary_line in gene_summary_file:
			items = summary_line.split()
			sample = items[0].strip()
			marker_coverage = float(items[5])
			marker_coverage_samples.append(sample)
			marker_coverages.append(marker_coverage)
		
		gene_summary_file.close()
		marker_coverage_samples = su.parse_merged_sample_names(marker_coverage_samples)
		marker_coverage_map = {sample: marker_coverage for sample,marker_coverage in zip(marker_coverage_samples, marker_coverages)}
		
		# Get rid of headers + get list of samples
		reads_line = gene_reads_file.readline() # header
		depth_line = gene_depth_file.readline() # header
		presabs_line = gene_presabs_file.readline() # header
		items = reads_line.split() # gene_id, samples...
		samples = su.parse_merged_sample_names(items[1:])
		
		# Restrict to intersection of allowed_samples with samples
		allowed_samples = set(samples) if (allowed_samples is None) else (set(allowed_samples) & set(samples))
		
		# Boolean array: whether each of allowed_samples is in samples
		desired_sample_idxs = numpy.array([(s in allowed_samples) for s in samples])
		
		# Marker coverages in same order as samples
		marker_coverages = numpy.array([marker_coverage_map[s] for s in samples])
		
		# Final version of samples and marker coverages
		desired_samples = samples[desired_sample_idxs]
		marker_coverages = marker_coverages[desired_sample_idxs]
		
		gene_presence_matrix = []
		gene_reads_matrix = []
		gene_depth_matrix = []
		gene_names = []
		
		num_genes_processed = 0
		
		# Read first lines
		reads_line = gene_reads_file.readline()
		depth_line = gene_depth_file.readline()
		presabs_line = gene_presabs_file.readline()
		
		while reads_line != "":
				
				# First check for gene presence/absence
				items = presabs_line.split()
				gene_name = items[0]
				gene_presences = numpy.array([float(item) for item in items[1:]])[desired_sample_idxs]
				
				if True: # gene_presences.sum() > 0.5:
				
						gene_reads = numpy.array([float(item) for item in reads_line.split()[1:]])[desired_sample_idxs]
						gene_depths = numpy.array([float(item) for item in depth_line.split()[1:]])[desired_sample_idxs]
						
						# Note to self: not uniform across samples!
						#gene_lengths = gene_reads/(gene_depths+(gene_reads<0.5))
						#print gene_lengths
						
						# gene is present in at least one individual! 
						gene_presence_matrix.append(gene_presences)
						gene_depth_matrix.append(gene_depths)
						gene_reads_matrix.append(gene_reads)
						gene_names.append(gene_name)		
				
				num_genes_processed+=1
				
				reads_line = gene_reads_file.readline() # header
				depth_line = gene_depth_file.readline() # header
				presabs_line = gene_presabs_file.readline() # header
		
		gene_reads_file.close()
		gene_depth_file.close()
		gene_presabs_file.close()
		gene_presence_matrix = numpy.array(gene_presence_matrix)
		gene_depth_matrix = numpy.array(gene_depth_matrix)
		gene_reads_matrix = numpy.array(gene_reads_matrix)

		if convert_centroid_names:
				new_gene_names = []
				centroid_gene_map = midas_db_utils.load_centroid_gene_map(species_name)
				for gene_name in gene_names:
						new_gene_names.append(centroid_gene_map[gene_name])
		else:
				new_gene_names=gene_names
		
		new_gene_names = numpy.array(new_gene_names)
				
		# Now weed out disallowed genes if provided
		disallowed_genes=set(disallowed_genes)
		allowed_gene_idxs = []
		for gene_idx in xrange(0,len(new_gene_names)):
				
				if new_gene_names[gene_idx] in disallowed_genes:
						# don't include
						pass
				else:
						allowed_gene_idxs.append(gene_idx)
		allowed_gene_idxs = numpy.array(allowed_gene_idxs)
		
		new_gene_names = new_gene_names[allowed_gene_idxs]
		gene_presence_matrix = gene_presence_matrix[allowed_gene_idxs,:]
		gene_depth_matrix = gene_depth_matrix[allowed_gene_idxs,:]
		gene_reads_matrix = gene_reads_matrix[allowed_gene_idxs,:]
		
		return desired_samples, new_gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix

