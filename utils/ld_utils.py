from utils import config
import gzip
import numpy
import os

intermediate_filename_template = '%s%s.txt.gz'
ld_directory = '%slinkage_disequilibria/' % (config.data_directory)

def load_ld_map(species_name):

    ld_map = {}

    intermediate_filename = intermediate_filename_template % (ld_directory, species_name)

    if not os.path.isfile(intermediate_filename):
        return ld_map

    file = gzip.open(intermediate_filename,"r")
    header_line = file.readline() # header
    header_items = header_line.split(",")
    
    distance_strs = [item.split(":")[-1] for item in header_items[4:]]
    
    distances = []
    intragene_idxs = []
    
    intergene_distances = []
    intergene_idxs = []
    
    control_idx = -1
    
    for i in xrange(0,len(distance_strs)-1):
        
        if distance_strs[i].startswith('g'):
            # an intergene distance
            intergene_idxs.append(i)
            intergene_distances.append(long(distance_strs[i][1:]))
        else:
            # an intragene distance
            intragene_idxs.append(i)
            distances.append(float(distance_strs[i]))
            
    distances = numpy.array(distances)
    intragene_idxs = numpy.array(intragene_idxs)
    
    intergene_distances = numpy.array(intergene_distances)
    intergene_idxs = numpy.array(intergene_idxs)
    
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
        
        clade_type = items[1].strip()
        variant_type = items[2].strip()
        pi = float(items[3])
        
        rsquared_numerators = []
        rsquared_denominators = []
        lds = []
        counts = []
        for item in items[4:]:
            subitems = item.split(":")
            rsquared_numerators.append(float(subitems[0]))
            rsquared_denominators.append(float(subitems[1]))
            counts.append(float(subitems[2]))
        
        rsquared_numerators = numpy.array(rsquared_numerators)
        rsquared_denominators = numpy.array(rsquared_denominators)
        counts = numpy.array(counts)
        
        
        lds = rsquared_numerators/rsquared_denominators
        
        control_numerator = rsquared_numerators[control_idx]
        control_denominator = rsquared_denominators[control_idx]
        control_count = counts[control_idx]
        
        control_ld = control_numerator/control_denominator
        
        intragene_rsquared_numerators = rsquared_numerators[intragene_idxs]
        intragene_rsquared_denominators = rsquared_denominators[intragene_idxs]
        intragene_counts = counts[intragene_idxs]
        
        intergene_rsquared_numerators = rsquared_numerators[intergene_idxs]
        intergene_rsquared_denominators = rsquared_denominators[intergene_idxs]
        intergene_counts = counts[intergene_idxs]
        
        ld_map[(clade_type, variant_type)] = (distances, intragene_rsquared_numerators, intragene_rsquared_denominators, intragene_counts, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_counts, control_numerator, control_denominator, control_count, pi)
        
    return ld_map