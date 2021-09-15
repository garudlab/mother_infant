import numpy, os.path, gzip
import config, parse_midas_data

# ===========================================================================
# This module (core_gene_utils) contains the following utilities:
# 
# 	parse_core_genes
# 	parse_shared_genes
# 	parse_non_shared_reference_genes
# 	parse_non_shared_pangenome_genes
# 	get_good_pangenome_samples
# 	parse_gene_freqs
#
# ===========================================================================

# Prevalence cohorts: adult, hmp, infant
# where all is union of core genes from the other three

def get_filename(type, prev_cohort, external=False, species=None):
	
	if external:
		core_genes_directory = "%s/core_genes/external/" % (config.data_directory)
	else:
		core_genes_directory = "%s/core_genes/" % (config.data_directory)
	
	prev_cohort_directory = '%s/%s/' % (core_genes_directory, prev_cohort)
	
	# Type must be one of shared_genes, core_genes, core_genes_stringent, or gene_freqs
	# default to core
	if type == 'gene_freqs' and species != None:
		return '%s/prev_%s/%s_%s.txt.gz' % (prev_cohort_directory, prev_cohort, species, type)
	elif type in ['shared_genes', 'core_genes', 'core_genes_stringent']:
		return '%s/%s.txt.gz' % (prev_cohort_directory, type)
	else:
		raise ValueError("Bad arguments for core genes filename")

# ===========================================================================
# Returns set of core genes for specified set of species
# ===========================================================================

def parse_core_genes(desired_species_name = None, prev_cohort='all', external_filtering=True):
		
		core_genes = set() # Core genes for the specified species
		
		if prev_cohort == 'all':
			# Take union of core genes from adult, hmp, infant
			for sub_prev_cohort in ['adult', 'hmp', 'infant']:
				core_gene_filename = get_filename('core_genes', sub_prev_cohort)
				core_gene_file = gzip.GzipFile(core_gene_filename,"r")
				for line in core_gene_file:
					items = line.split(':')
					if len(items)<2:
						continue
					species_name = items[0].strip()
					gene_names = [subitem.strip() for subitem in items[1].split(",")]		
					if (species_name == desired_species_name) or (desired_species_name == None):
						core_genes.update(gene_names)
				
				core_gene_file.close()
			
		else:
			# Otherwise only consider specified prevalence cohort
			core_gene_filename = get_filename('core_genes', prev_cohort)
			core_gene_file = gzip.GzipFile(core_gene_filename,"r")
			for line in core_gene_file:					
					items = line.split(":")
					if len(items)<2:
							continue
					species_name = items[0].strip()
					gene_names = [subitem.strip() for subitem in items[1].split(",")]					
					if (species_name == desired_species_name) or (desired_species_name == None):
							core_genes.update(gene_names)
					
			core_gene_file.close()
		
		# Account for externally provided core genes if available
		external_core_gene_filename = get_filename('core_genes', prev_cohort, external=True)
		external_core_genes = set()
		
		if os.path.isfile(external_core_gene_filename):				
				external_core_gene_file = gzip.GzipFile(external_core_gene_filename,"r")				
				for line in external_core_gene_file:				
						items = line.split(":")
						if len(items)<2:
								continue						
						species_name = items[0].strip()
						gene_names = [subitem.strip() for subitem in items[1].split(",")]				
						if (species_name == desired_species_name) or (desired_species_name == None):
								external_core_genes.update(gene_names)						
				external_core_gene_file.close() 
		
		if external_filtering and len(external_core_genes)>0:
				core_genes = (core_genes & external_core_genes)
		
		return core_genes

# ===========================================================================
# Returns set of shared genes for specified set of species
# ===========================================================================

def parse_shared_genes(desired_species_name = None, prev_cohort='all', external_filtering = True):
		
		shared_genes = set()
		
		if prev_cohort == 'all':
			
			# Take union of core genes from adult, hmp, infant
			for sub_prev_cohort in ['adult', 'hmp', 'infant']:
				
				shared_gene_filename = get_filename('shared_genes', sub_prev_cohort)
				shared_gene_file = gzip.GzipFile(shared_gene_filename,"r")
				
				for line in shared_gene_file:
					items = line.split(":")
					
					if len(items)<2:
						continue
					
					species_name = items[0].strip()
					gene_names_str = items[1].strip()					
					# N/A means wasn't enough pangenome data to detect shared genes
					gene_names = [] if gene_names_str.startswith('N/A') else [subitem.strip() for subitem in gene_names_str.split(",")]
					
					if (species_name==desired_species_name) or (desired_species_name==""):
							shared_genes.update(gene_names)
				
				shared_gene_file.close() 
			
		else:
			# Otherwise only consider specified prevalence cohort
			shared_gene_filename = get_filename('shared_genes', prev_cohort)
			shared_gene_file = gzip.GzipFile(shared_gene_filename,"r")
			
			for line in shared_gene_file:
				items = line.split(":")
				
				if len(items)<2:
					continue
				
				species_name = items[0].strip()
				gene_names_str = items[1].strip()					
				# N/A means wasn't enough pangenome data to detect shared genes
				gene_names = [] if gene_names_str.startswith('N/A') else [subitem.strip() for subitem in gene_names_str.split(",")]
				
				if (species_name==desired_species_name) or (desired_species_name==""):
						shared_genes.update(gene_names)
			
			shared_gene_file.close()
		
		external_shared_gene_filename = get_filename('shared_genes', prev_cohort, external=True)
		external_shared_genes = set()
		
		if os.path.isfile(external_shared_gene_filename):
				
				external_shared_gene_file = gzip.GzipFile(external_shared_gene_filename,"r")
				
				for line in external_shared_gene_file:
						
						items = line.split(":")
						if len(items)<2:
								continue
						
						species_name = items[0].strip()
						gene_names_str = items[1].strip()
						if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
								gene_names = []
						else:
								gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
				
						if (species_name==desired_species_name) or (desired_species_name==""):
								external_shared_genes.update(gene_names)
						
				external_shared_gene_file.close() 
		
		if external_filtering and len(external_shared_genes)>0:
				# some externally provided core genes
				shared_genes = (shared_genes | external_shared_genes)
				
		return shared_genes

# ===========================================================================
# Returns set of reference genes which are not shared
# ===========================================================================

def parse_non_shared_reference_genes(desired_species_name="", prev_cohort='all', external_filtering=True):
		
		from utils import parse_midas_data
		
		shared_genes = parse_shared_genes(desired_species_name, prev_cohort, external_filtering)
		reference_genes = parse_midas_data.load_reference_genes(desired_species_name)
		non_shared_reference_genes = set(reference_genes)-shared_genes
		
		return non_shared_reference_genes

# ===========================================================================
# Returns set of pangenome genes which are not shared
# ===========================================================================

def parse_non_shared_pangenome_genes(desired_species_name="", prev_cohort='all', external_filtering=True):
		
		from utils import parse_midas_data
		
		shared_genes = parse_shared_genes(desired_species_name, prev_cohort, external_filtering)
		pangenome_genes, pangenome_centroid_genes = parse_midas_data.load_pangenome_genes(desired_species_name)
		# TODO: Not sure if I should be using the first or second
		non_shared_pangenome_genes = set(pangenome_centroid_genes)-shared_genes
		return non_shared_pangenome_genes

# ===========================================================================
# Returns indices for samples which have enough present genes (copy number
# exceeds a low threshold) of which not too many are high-copynum
# ===========================================================================

def get_good_pangenome_samples(marker_coverages, gene_copynum_matrix, species_name):
		
		cmin = config.core_genome_min_copynum
		cmax = config.core_genome_max_copynum	 
		
		# For each sample, get number of genes which are "present" (copynum > 0.3)
		# and which are "high" (copynum > 3)
		num_present_genes = (gene_copynum_matrix>cmin).sum(axis=0)
		num_high_genes = (gene_copynum_matrix>cmax).sum(axis=0)
		
		# Get number of reference genes
		# TODO: does not match number of genes in copynum matrix
		num_reference_genes = len(parse_midas_data.load_reference_genes(species_name))
		
		# Want at least 30% of all reference genes to be present, and want no more
		# than 30% of present genes to be high
		min_present_genes = 0.3*num_reference_genes
		max_high_genes = 0.3*num_present_genes
		
		good_sample_idxs = (num_present_genes>min_present_genes)*(num_high_genes<max_high_genes)
		
		return good_sample_idxs

# ===========================================================================
# Returns gene frequency map (gene name -> prevalence)
# ===========================================================================

def parse_gene_freqs(desired_species_name, prev_cohort='hmp', use_external=False):
		
		filename = get_filename('gene_freqs', prev_cohort, external=use_external, species=desired_species_name)
		
		if not os.path.isfile(filename):
				return {}
				
		file = gzip.open(filename,"r")
		gene_freq_map = {}
		for line in file:
				items = line.split()
				gene_name = items[0]
				f = float(items[1])
				gene_freq_map[gene_name] = f
		file.close()
		
		return gene_freq_map
