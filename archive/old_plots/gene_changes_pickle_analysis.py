import pickle
from collections import defaultdict
from operator import itemgetter

import config
import parse_midas_data

cohort = 'hmp' # hmp or infant

good_species_list = parse_midas_data.parse_good_species_list()

# Dictionaries for storing pickled data
gene_snp_change_dict = defaultdict(int)
gene_gainloss_dict = defaultdict(int)
gene_change_species = defaultdict(list)

# Summary output files
output_fn_sp = open("%s/gene_changes_summary_by_species.txt" % config.data_directory, "w")
output_fn_sum = open("%s/gene_changes_summary.txt" % config.data_directory, "w")

for species in good_species_list:
	try:
		gene_changes_pkl = pickle.load(open("%s/pickles/gene_changes_%s/%s_gene_changes.p" % (config.data_directory, cohort, species), 'r'))
	except:
		print("Couldn't open gene changes pickle for %s" % species)
		continue
	
	output_fn_sp.write("%s\tsnps\tall\n" % species)
	
	gene_changes = gene_changes_pkl[species]
	if gene_changes != {}:
		gene_changes = gene_changes['gene_changes']
	else:
		continue
	
	for gene_name in gene_changes.keys():
		num_snp_changes = gene_changes[gene_name]['snps']
		gene_snp_change_dict[gene_name] += num_snp_changes
		
		num_gene_gainloss = gene_changes[gene_name]['all']
		gene_gainloss_dict[gene_name] += num_gene_gainloss
		
		gene_change_species[gene_name].append(species)
		
		output_fn_sp.write(gene_name)
		output_fn_sp.write('\t')
		output_fn_sp.write(str(num_snp_changes))
		output_fn_sp.write('\t')
		output_fn_sp.write(str(num_gene_gainloss))
		output_fn_sp.write('\n')
	
	output_fn_sp.write('\n')

output_fn_sum.write('Gene name\tSnp changes\tGene gain/loss\nNumber Species\n')
sorted_gene_changes = sorted(gene_snp_change_dict.items(), reverse=True, key=itemgetter(1))

for gene_name, num_snp_changes in sorted_gene_changes:
	output_fn_sum.write(gene_name)
	output_fn_sum.write('\t')
	output_fn_sum.write(str(num_snp_changes))
	output_fn_sum.write('\t')
	output_fn_sum.write(str(gene_gainloss_dict[gene_name]))
	output_fn_sum.write('\t')
	output_fn_sum.write(str(len(gene_change_species[gene_name])))
	output_fn_sum.write('\n')
