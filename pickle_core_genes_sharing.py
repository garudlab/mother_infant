from utils import parse_midas_data as pmd, core_gene_utils as cgu
import config

good_species_list = pmd.parse_good_species_list()

f = open('%s/core_genes_sharing.csv' % config.analysis_directory, 'w')

f.write((',').join(['species', 'num_total', 'num_hmp', 'num_hmp_infant_shared', 'num_infant', 'num_adult_infant_shared', 'num_adult']))
f.write('\n')

for species in good_species_list:
	print('Working on %s...' % species)
	core_genes_all = cgu.parse_core_genes(species, prev_cohort='all')
	core_genes_hmp = cgu.parse_core_genes(species, prev_cohort='hmp')
	core_genes_infant = cgu.parse_core_genes(species, prev_cohort='infant')
	core_genes_adult = cgu.parse_core_genes(species, prev_cohort='adult')
	
	core_genes_hmp_infant_shared = core_genes_hmp.intersection(core_genes_infant)
	core_genes_adult_infant_shared = core_genes_adult.intersection(core_genes_infant)
	
	data_items = [species, len(core_genes_all), len(core_genes_hmp), len(core_genes_hmp_infant_shared), len(core_genes_infant), len(core_genes_adult_infant_shared), len(core_genes_adult)]
	
	row_str = (',').join([str(item) for item in data_items])
	f.write(row_str + '\n')

f.close()
