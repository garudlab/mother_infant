import matplotlib	 
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import parse_midas_data as pmd
import midas_db_utils
from collections import defaultdict

def get_species_color_map(all_species = pmd.parse_species_list()):
	
	gfo_phylum_map = midas_db_utils.load_gfo_phylum_map()
	
	all_phyla = set()
	phylum_species_map = defaultdict(list)
	
	for species in all_species:
		gfo = species.split('_')[0] # Either genus, family or order
		try:
			phylum = gfo_phylum_map[gfo]
			all_phyla.add(phylum)
			phylum_species_map[phylum].append(species)
		except:
			print(species) # Guyana massiliensis is unclassified, skip for now
	
	# There are the following 7 phyla in these gut microbiota
	species_color_map = {}
	ordered_species_list = []
	
	# Firmicutes (81)
	species_list = sorted(phylum_species_map['Firmicutes'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.PuRd, len(species_list), 7, 3)	
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)
	
	# Bacteroidetes (54)
	species_list = sorted(phylum_species_map['Bacteroidetes'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.YlGnBu, len(species_list), 5, 2)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)
	
	# Proteobacteria (10)
	species_list = sorted(phylum_species_map['Proteobacteria'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.Oranges, len(species_list), 2, 3)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)
	
	# Actinobacteria (9)
	species_list = sorted(phylum_species_map['Actinobacteria'], key=lambda x: x[1])
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.Wistia, len(species_list), 2, 4)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)
	
	# Fusobacteria (1), Spirochaetes (1), Verrucomicrobia (1)
	species_list = phylum_species_map['Fusobacteria'] + phylum_species_map['Spirochaetes'] + phylum_species_map['Verrucomicrobia']
	ordered_species_list += species_list
	m = get_cm_ScalerMappable(matplotlib.cm.Greens, len(species_list), 2, 1)
	for i in range(len(species_list)):
		species_color_map[species_list[i]] = m.to_rgba(i)
	
	# Special: set '-' to gray
	species_color_map['-'] = (0.6, 0.6, 0.6, 1.0)
	
	return species_color_map, ordered_species_list

def get_cm_ScalerMappable(cmap, num_colors, offset1 = 0, offset2 = 0):
	norm = matplotlib.colors.Normalize(vmin = 0 - offset1, vmax = num_colors- 1 + offset2)
	return matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)

def list_to_colors(input_list):
	cmap = matplotlib.cm.hsv
	
	input_list_dict = {}
	i = 0
	for elem in set(input_list):
		input_list_dict[elem] = i
		i += 1
	
	norm = matplotlib.colors.Normalize(vmin=0, vmax=i-1)
	m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
	
	color_list = []
	for elem in input_list:
		color_i = input_list_dict[elem]
		color_list.append(m.to_rgba(color_i))
	
	return color_list

def colors_to_legend_elements(colors, labels):
	legend_elements = {}
	for color, label in zip(colors, labels):
		legend_elements[color] = matplotlib.patches.Patch(facecolor=color, label=label)
	return [legend_elements[c] for c in remove_list_duplicates(colors)]

def remove_list_duplicates(mylist):
	# Returns list with only unique elements (like set)
	# but ordered according to first appearnce in original list
	result = []
	for elem in mylist:
		if elem not in result:
			result.append(elem)
	return result


def get_pretty_species_name(species_name, include_number=False):
	items = species_name.split("_")	
	pretty_name = "%s %s" % (items[0], items[1])	
	if include_number:
			pretty_name += (" (%s)" % (items[2]))			
	return pretty_name
    
def get_abbreviated_species_name(species_name):    
	items = species_name.split("_")	
	pretty_name = "%s. %s" % (items[0][0], items[1])			
	return pretty_name
