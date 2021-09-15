import numpy
import gzip
import os
import config

substitution_rate_directory = '%s/substitution_rates/' % (config.data_directory)
intermediate_filename_template = '%s/%s/%s.txt.gz'

def load_substitution_rate_map(species_name, prev_cohort='all'):

    intermediate_filename = intermediate_filename_template % (substitution_rate_directory, prev_cohort, species_name)

    substitution_rate_map = {}

    if not os.path.isfile(intermediate_filename):
        return substitution_rate_map
    
    file = gzip.open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        record_strs = [", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'Num_muts', 'Num_revs', 'Num_mut_opportunities', 'Num_rev_opportunities'])]    
            
        sample_1 = items[1].strip()
        sample_2 = items[2].strip()
        type = items[3].strip()
        num_muts = float(items[4])
        num_revs = float(items[5])
        num_mut_opportunities = float(items[6])
        num_rev_opportunities = float(items[7])
        
        num_changes = num_muts+num_revs
        num_opportunities = num_mut_opportunities+num_rev_opportunities
        
        sample_pair = (sample_1, sample_2)
        
        if type not in substitution_rate_map:
            substitution_rate_map[type] = {}
          
        substitution_rate_map[type][sample_pair] = (num_muts, num_revs, num_mut_opportunities, num_rev_opportunities)
        
    return substitution_rate_map

def calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, type, allowed_samples=[]): 
    # Rewritten to preserve order of allowed samples
    # If allowed samples contains things that are not in DB, it returns zero opportunities

    total_sample_set = set([])
    for sample_1, sample_2 in substitution_rate_map[type].keys():
        total_sample_set.add(sample_1)
        total_sample_set.add(sample_2)

    if len(allowed_samples)==0:
        allowed_samples = list(sorted(total_sample_set))    
        
    # allows us to go from sample name to idx in allowed samples (to preserve order)
    sample_idx_map = {allowed_samples[i]:i for i in xrange(0,len(allowed_samples))}
    
    mut_difference_matrix = numpy.zeros((len(allowed_samples), len(allowed_samples)))*1.0
    rev_difference_matrix = numpy.zeros_like(mut_difference_matrix)
    
    mut_opportunity_matrix = numpy.zeros_like(mut_difference_matrix)
    rev_opportunity_matrix = numpy.zeros_like(mut_difference_matrix)
    
    for sample_pair in substitution_rate_map[type].keys():
        
        sample_i = sample_pair[0]
        sample_j = sample_pair[1]
        
        if not ((sample_i in sample_idx_map) and (sample_j in sample_idx_map)):
            continue
        
        i = sample_idx_map[sample_i]
        j = sample_idx_map[sample_j]
        
        num_muts, num_revs, num_mut_opportunities, num_rev_opportunities = substitution_rate_map[type][sample_pair]
        
        mut_difference_matrix[i,j] = num_muts
        rev_difference_matrix[i,j] = num_revs
        
        mut_opportunity_matrix[i,j] = num_mut_opportunities
        rev_opportunity_matrix[i,j] = num_rev_opportunities
        
    return allowed_samples, mut_difference_matrix, rev_difference_matrix, mut_opportunity_matrix, rev_opportunity_matrix

    
def calculate_matrices_from_substitution_rate_map(substitution_rate_map, type, allowed_samples=[]):
# once the map is loaded, then we can compute rate matrices in this definition (so, it relies on the previous def)    

    samples, mut_difference_matrix, rev_difference_matrix, mut_opportunity_matrix, rev_opportunity_matrix = calculate_mutrev_matrices_from_substitution_rate_map( substitution_rate_map, type, allowed_samples)

    difference_matrix = mut_difference_matrix+rev_difference_matrix
    opportunity_matrix = mut_opportunity_matrix+rev_opportunity_matrix
    
    return samples, difference_matrix, opportunity_matrix
