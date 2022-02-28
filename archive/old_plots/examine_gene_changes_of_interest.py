import sample_utils as su, config, parse_midas_data
import diversity_utils
import pickle
import sys

def format_tp_label(cohort, tp_pair):
	tpa, tpb = tp_pair
	if tpa[0] == 'M' and tpb[0] == 'I' or tpa[0] == 'I' and tpb[0] == 'M':
		tp1, tp2 = (tpa, tpb) if tpa[0] == 'M' else (tpb, tpa)
	else:
		tp1, tp2 = (tpa, tpb) if (tpa < tpb) else (tpb, tpa)
	subj1, subj2 = tp1[0], tp2[0]
	order1, order2 = int(tp1[1]), int(tp2[1])
	if cohort == 'backhed':
		ilabels = ['I_Birth', 'I_4mon','I_12mon']
		mlabels = ['M_Birth']
	elif cohort == 'ferretti':
		ilabels = ['I_1day','I_3days', 'I_1wk', 'I_1mon', 'I_4mon']
		mlabels = ['M_Birth']
	elif cohort == 'yassour':
		ilabels = ['I_Birth', 'I_2wk','I_1mon', 'I_2mon', 'I_3mon']
		mlabels = ['M_Gest', 'M_Birth', 'M_3mon']
	elif cohort == 'hmp':
		return (str(order1) + " > " + str(order2))
	else:
		return 'error: bad cohort'
	label1 = mlabels[order1-1] if subj1 == 'M' else ilabels[order1-1]
	label2 = mlabels[order2-1] if subj2 == 'M' else ilabels[order2-1]
	return (label1.replace('_','-') + " > " + label2.replace('_','-'))

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = su.parse_subject_sample_map()
sample_order_map = su.parse_sample_order_map()
sample_country_map = su.parse_sample_country_map()
sample_subject_map = su.calculate_sample_subject_map(subject_sample_map)
sys.stderr.write("Done!\n")

all_samples = sample_country_map.keys()

# HMP samples
hmp_samples = [x for x in all_samples if sample_country_map[x] == 'United States']

# Mother-infant samples
mother_samples = su.get_sample_names('mother','all')
infant_samples = su.get_sample_names('infant','all')
backhed_samples = su.get_sample_names('backhed', 'all')
ferretti_samples = su.get_sample_names('ferretti', 'all')
yassour_samples = su.get_sample_names('yassour', 'all')
	
# Investigate genes of interest
genes_of_interest = pickle.load(open("%s/pickles/gene_changes_hmp/gene_changes_of_interest.pkl" % config.data_directory, "r"))

# Genome ID -> species map
genome_species_dict = {}
good_species_list = parse_midas_data.parse_good_species_list()
for species in good_species_list:
	for genome_id in parse_midas_data.get_ref_genome_ids(species):
		genome_species_dict[genome_id] = species

for gene_name in genes_of_interest:
	print gene_name
	
	subjects = []
	species = []
	for sample_i, sample_j in genes_of_interest[gene_name]:
		subject = sample_order_map[sample_i][0][:-2]
		subjects.append(subject)
		
		if sample_i and sample_j in hmp_samples:
			tp_i, tp_j = ('A' + str(sample_order_map[sample_i][1]), 'A' + str(sample_order_map[sample_j][1]))
			cohort = 'hmp'
		else:
			tp_i = ('M' if sample_i in mother_samples else 'I') + str(sample_order_map[sample_i][1])
			tp_j = ('M' if sample_j in mother_samples else 'I') + str(sample_order_map[sample_j][1])
			if sample_i and sample_j in backhed_samples:
				cohort = 'backhed'
			elif sample_i and sample_j in ferretti_samples:
				cohort = 'ferretti'
			elif sample_i and sample_j in yassour_samples:
				cohort = 'yassour'
		
		tp_pair = frozenset((tp_i, tp_j))
		
		gene_changes = genes_of_interest[gene_name][frozenset((sample_i, sample_j))]
		for i in range(len(gene_changes)):
			print format_tp_label(cohort, tp_pair) # Duplicate print for multiple gene changes in the same sample pair
			gene_change = gene_changes[i]
			gene_id = gene_change[0]
			if len(gene_change) == 5: # Gene gain/loss: (gene_name, D1, Dm1, D2, Dm2)
				_, D1, Dm1, D2, Dm2 = gene_change
				print("Coverage: %i %i %i %i" % (D1, Dm1, D2, Dm2))
			elif len(gene_change) == 8: # SNP change: (gene_name, contig, position, variant_type, A1, D1, A2, D2)
				_, contig, position, variant_type, A1, D1, A2, D2 = gene_change
				print("Variant type: " + variant_type)
				print("Coverage: %i %i %i %i" % (A1, D1, A2, D2))
		
		genome_id = gene_id.split('.')[0] + "." + gene_id.split('.')[1]
		species_name = genome_species_dict[genome_id]
		species.append(species_name)
		
		# dummy_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D']))
		
		# perr = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j)[0]
		
		# print("Error rate: " + str(perr))
		print("------------------------------------------------")
	
	print("All species: " + str(species))
	print("All subjects: " + str(subjects) + "\n")	
