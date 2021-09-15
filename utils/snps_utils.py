import config
import numpy
import gzip
import os

def parse_preexisting_snps(species_name):
	
	intermediate_filename = '%s/preexisting_snps.txt.gz' % (config.data_directory)		
	
	preexisting_snps = {}	 
	file = gzip.GzipFile(intermediate_filename,"r")
	for line in file:
			if line.startswith(species_name):
					contig_items = line.split(";")[1:]
					for contig_item in contig_items:
							contig_item = contig_item.strip()
							if contig_item=="":
									continue
							contig_subitems = contig_item.split(":")
							contig = contig_subitems[0].strip()
							snp_items = contig_subitems[1].split()
							for snp_item in snp_items:
									snp_subitems = snp_item.split(",")
									location = long(snp_subitems[0])
									prevalence = float(snp_subitems[1])
									
									preexisting_snps[(contig,location)] = prevalence
									
	file.close()
	
	return preexisting_snps

def load_private_snv_map(species_name, prev_cohort='all'):
	
	private_snv_directory = '%s/private_snvs/%s/' % (config.data_directory, prev_cohort)
	intermediate_filename_template = '%s/%s.txt.gz'	
	intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)

	private_snv_map = {}

	if not os.path.isfile(intermediate_filename):
			return private_snv_map
	
	file = gzip.open(intermediate_filename,"r")
	file.readline() # header
	for line in file:
	 
			items = line.split(",")
			
			contig = items[0].strip()
			location = long(items[1])
			gene_name = items[2].strip()
			variant_type = items[3].strip()
			host = items[4].strip()
			
			private_snv_map[(contig, location)] = (gene_name, variant_type, host)
			
	return private_snv_map

def load_snv_distance_map(species_name):
		
		private_snv_directory = '%s/snv_distances/' % (config.data_directory)
		intermediate_filename_template = '%s/%s.txt.gz'	
		intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)

		snv_distance_map = {}

		file = gzip.open(intermediate_filename,"r")
		file.readline() # header
		for line in file:
		 
				items = line.split(",")
				
				contig = items[0].strip()
				location = long(items[1])
				variant_type = items[2].strip()
				derived_allele_count = long(items[3])
				ancestral_allele_count = long(items[4])
				min_between_d = float(items[5])
				max_within_d1 = float(items[6])
				max_within_d2 = float(items[7])
				
				snv_distance_map[(contig, location)] = (variant_type, derived_allele_count, ancestral_allele_count, min_between_d, max_within_d1, max_within_d2)
				
		return snv_distance_map

# Loading file
def parse_snp_prevalences(cohort, desired_species_name):
		
		intermediate_filename_template = config.data_directory+"snp_prevalences_%s/%s.txt.gz"
		intermediate_filename_alt_template = config.data_directory+"snp_prevalences_alt_%s/%s.txt.gz"
		
		if cohort == 'all':
			intermediate_filename = intermediate_filename_template % (cohort, desired_species_name)
		else:
			intermediate_filename = intermediate_filename_alt_template % (cohort, desired_species_name)
				
		snp_prevalences = {}
		
		if not os.path.isfile(intermediate_filename):
				return snp_prevalences
		
		file = gzip.GzipFile(intermediate_filename,"r")
		file.readline()
		for line in file:
				items = line.split(",")
				contig = items[0]
				location = long(items[1])
				population_freq = float(items[2])
				snp_freq = float(items[3])
				
				snp_prevalences[(contig,location)] = snp_freq
														
		file.close()
		
		return snp_prevalences
		
# originally from calculate_snp_prevalences
def parse_population_freqs(cohort, desired_species_name, polarize_by_consensus=False):
		
		intermediate_filename_template = config.data_directory+"snp_prevalences/%s/%s.txt.gz"
		intermediate_filename_alt_template = config.data_directory+"snp_prevalences_alt/%s/%s.txt.gz"
		
		if cohort == 'all':
			intermediate_filename = intermediate_filename_template % (cohort, desired_species_name)
		else:
			intermediate_filename = intermediate_filename_template % (cohort, desired_species_name)
			
		population_freqs = {}
		
		if not os.path.isfile(intermediate_filename):
				return population_freqs

		file = gzip.GzipFile(intermediate_filename,"r")
		file.readline()
		for line in file:
				items = line.split(",")
				contig = items[0]
				location = long(items[1])
				population_freq = float(items[2])
				snp_freq = float(items[3])
				
				if polarize_by_consensus:
						if population_freq > 0.5:
								population_freq = 1-population_freq
				
				if population_freq==0:
						pass
				else:
						population_freqs[(contig,location)] = population_freq
		
		file.close()
		
		return population_freqs

def load_singleton_rate_map(species_name):
# This definition is called whenever another script downstream uses the output of this data.
		
		singleton_directory = '%s/singleton_rates/' % (config.data_directory)
		intermediate_filename_template = '%s/%s.txt.gz'		
		intermediate_filename = intermediate_filename_template % (singleton_directory, species_name)

		singleton_rate_map = {}

		if not os.path.isfile(intermediate_filename):
				return singleton_rate_map
		
		file = gzip.open(intermediate_filename,"r")
		file.readline() # header
		for line in file:
				items = line.split(",")
				if items[0].strip()!=species_name:
						continue
						
				sample_i = items[1].strip()
				sample_j = items[2].strip()
				type = items[3].strip()
				num_singletons = float(items[4])
				num_doubletons = float(items[5])
				num_differences = float(items[6])
				num_opportunities = float(items[7])
				
				if type not in singleton_rate_map:
						singleton_rate_map[type] = {}
				
				if sample_i==sample_j:
						num_singletons = 0
						num_doubletons = 0
						num_differences = 0
					
				singleton_rate_map[type][sample_i, sample_j] = (num_singletons, num_doubletons, num_differences, num_opportunities)
				
		return singleton_rate_map

def calculate_matrices_from_singleton_rate_map(singleton_rate_map, type, allowed_samples=[]):
# once the map is loaded, then we can compute rate matrices in this definition (so, it relies on the previous def)		

		sample_set = set([])
		for sample in singleton_rate_map[type].keys():
				sample_set.add(sample)
		
		if len(allowed_samples)>0:
				allowed_sample_set = set(allowed_samples)		 
		else:
				allowed_sample_set = sample_set
		
		sample_set = set()
		for sample_i, sample_j in singleton_rate_map[type]:
				sample_set.add(sample_i)
				sample_set.add(sample_j)		
				 
		if len(allowed_samples)==0:
				allowed_samples = list(sorted(allowed_sample_set))
		
		
		samples = []
		# preserve same order as allowed samples
		for sample in allowed_samples:
				if sample in sample_set:
						samples.append(sample)
							 
		singleton_matrix = numpy.zeros((len(samples),len(samples)))*1.0
		doubleton_matrix = numpy.zeros_like(singleton_matrix)
		difference_matrix = numpy.zeros_like(singleton_matrix)
		opportunity_matrix = numpy.zeros_like(singleton_matrix)
		
		for i in xrange(0,len(samples)):
				for j in xrange(0,len(samples)):
						
						num_singletons, num_doubletons, num_differences, num_opportunities = singleton_rate_map[type][(samples[i], samples[j])]
				
						if i==j:
								num_doubletons = 0
				
						singleton_matrix[i,j] = num_singletons
						doubleton_matrix[i,j] = num_doubletons
						difference_matrix[i,j] = num_differences
						opportunity_matrix[i,j] = num_opportunities
		 
		#print singleton_matrix, opportunity_matrix		 
		return samples, singleton_matrix, doubleton_matrix, difference_matrix, opportunity_matrix

		
def calculate_matrices_from_substitution_rate_map(substitution_rate_map, type, allowed_samples=[]):
# once the map is loaded, then we can compute rate matrices in this definition (so, it relies on the previous def)		

		samples, mut_difference_matrix, rev_difference_matrix, mut_opportunity_matrix, rev_opportunity_matrix = calculate_mutrev_matrices_from_substitution_rate_map( substitution_rate_map, type, allowed_samples)

		difference_matrix = mut_difference_matrix+rev_difference_matrix
		opportunity_matrix = mut_opportunity_matrix+rev_opportunity_matrix
		
		return samples, difference_matrix, opportunity_matrix
