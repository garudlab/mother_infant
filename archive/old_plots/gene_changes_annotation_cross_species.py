# Within-host gene change annotation

import matplotlib as mpl
mpl.use('Agg')
import sample_utils as su, config, parse_midas_data, sample_utils
import pylab, sys, numpy, random
import diversity_utils, gene_diversity_utils, sfs_utils, calculate_substitution_rates, calculate_temporal_changes, parse_patric, species_phylogeny_utils, stats_utils
from collections import defaultdict
import operator
import numpy as np
import pickle

#######################################
# code for loading cross-species data #
#######################################

cohort = 'infant'
num_trials = 100

good_species_list = parse_midas_data.parse_good_species_list()

all_data = {} # Includes changes and nulls
gene_snp_change_dict = defaultdict(int) # Keys are genes
gene_gainloss_dict = defaultdict(int) # Keys are genes
gene_all_change_dict = defaultdict(int) # Keys are genes

gene_changes_file = open("%s/pickles/gene_changes_%s/all_species_gene_changes.p" % (config.data_directory, cohort), 'rb')
species_data = pickle.load(gene_changes_file)

for species_name in good_species_list:
	'''
	try:
		gene_changes_file = open("%s/pickles/gene_changes_%s/%s_gene_changes.p" % (config.data_directory, cohort, species_name), 'rb')
		species_data = pickle.load(gene_changes_file)
		if (len(species_data[species_name]) > 0):
			all_data[species_name] = species_data[species_name]
	except:
		sys.stderr.write("Missing gene changes file for " + species_name + '\n')
		continue
	'''
	if species_name in species_data:
		gene_changes = species_data[species_name]
	else:
		continue
	
	if gene_changes != {}:
		all_data[species_name] = gene_changes
		gene_changes = gene_changes['gene_changes']		
	else:
		continue
	
	for gene_name in gene_changes.keys():
		num_snp_changes = gene_changes[gene_name]['snps']
		if num_snp_changes != 0:
			gene_snp_change_dict[gene_name] += num_snp_changes
		
		num_gene_gainloss = gene_changes[gene_name]['all']
		if num_gene_gainloss != 0:
			gene_gainloss_dict[gene_name] += num_gene_gainloss
		
		gene_all_change_dict[gene_name] += (num_gene_gainloss + num_snp_changes)

# First look at individual genes: sort by number of changes
sorted_gene_all_change = sorted(gene_all_change_dict.items(), reverse=True, key=operator.itemgetter(1))

for gene_name, change_count in sorted_gene_all_change:
	if change_count > 1:
		print(str(change_count) + " " + gene_name)

print("Total number of gene changes: " + str(sum(gene_all_change_dict.values())))

sorted_gene_snp_change = sorted(gene_snp_change_dict.items(), reverse=True, key=operator.itemgetter(1))

for gene_name, change_count in sorted_gene_snp_change:
	if change_count > 1:
		between_host_null_exp = 0
		within_host_null_exp = 0
		for species in all_data:
			between_host_null = all_data[species]['null']['between_host_genes']
			within_host_null = all_data[species]['null']['present_genes']		
			if gene_name in between_host_null:
				between_host_null_exp += np.sum(between_host_null[gene_name]['snps'])/float(len(between_host_null[gene_name]['snps']))
			if gene_name in within_host_null:
				within_host_null_exp += np.sum(within_host_null[gene_name]['snps'])/float(len(within_host_null[gene_name]['snps']))
		print("%i	| %.2f | %.2f | %s" % (change_count, between_host_null_exp, within_host_null_exp, gene_name))

print("Total number of gene SNP changes: " + str(sum(gene_snp_change_dict.values())))

sorted_gene_gainloss = sorted(gene_gainloss_dict.items(), reverse=True, key=operator.itemgetter(1))

for gene_name, change_count in sorted_gene_gainloss:
	if change_count > 1:
		between_host_null_exp = 0
		within_host_null_exp = 0
		for species in all_data:
			between_host_null = all_data[species]['null']['between_host_genes']
			within_host_null = all_data[species]['null']['present_genes']		
			if gene_name in between_host_null:
				between_host_null_exp += np.sum(between_host_null[gene_name]['all'])/float(len(between_host_null[gene_name]['all']))
			if gene_name in within_host_null:
				within_host_null_exp += np.sum(within_host_null[gene_name]['all'])/float(len(within_host_null[gene_name]['all']))
		print("%i	| %.2f | %.2f | %s" % (change_count, between_host_null_exp, within_host_null_exp, gene_name))

print("Total number of gene gains/losses: " + str(sum(gene_gainloss_dict.values())))

# Get error information

perr_by_gene_snp_pkl_fn = '%s/pickles/perr_by_gene_snp_change.pkl' % config.data_directory
nerr_by_gene_snp_pkl_fn = '%s/pickles/nerr_by_gene_snp_change.pkl' % config.data_directory
perr_by_gene_gainloss_pkl_fn = '%s/pickles/perr_by_gene_gainloss.pkl' % config.data_directory
nerr_by_gene_gainloss_pkl_fn = '%s/pickles/nerr_by_gene_gainloss.pkl' % config.data_directory

perr_all_snp_changes = pickle.load(open(perr_by_gene_snp_pkl_fn, 'r'))
nerr_all_snp_changes = pickle.load(open(nerr_by_gene_snp_pkl_fn, 'r'))
perr_all_gene_gainloss = pickle.load(open(perr_by_gene_gainloss_pkl_fn, 'r'))
nerr_all_gene_gainloss = pickle.load(open(nerr_by_gene_gainloss_pkl_fn, 'r'))

gene_snp_change_errs = []
for desired_gene_name in gene_snp_change_dict:
	for gene_id in perr_all_snp_changes[cohort]:
		if parse_patric.load_patric_gene_name(gene_id) == desired_gene_name:
			print("yea")
			gene_snp_change_errs.append((desired_gene_name, perr_all_snp_changes[cohort][gene_id], nerr_all_snp_changes[cohort][gene_id]))

gene_gainloss_errs = []
for desired_gene_name in gene_gainloss_dict:
	for gene_id in perr_all_snp_changes[cohort]:
		if parse_patric.load_patric_gene_name(gene_id) == desired_gene_name:
			print("yea")
			gene_gainloss_errs.append((desired_gene_name, perr_all_gene_gainloss[cohort][gene_id], nerr_all_gene_gainloss[cohort][gene_id]))

# Sum gene changes across species
change_types = ['snps', 'gains', 'losses', 'all']
all_data['all_species'] = {'gene_changes': {}}

for species_name in all_data:
	if species_name != 'all_species':
		for gene in all_data[species_name]['gene_changes']:
			if gene not in all_data['all_species']['gene_changes']:
				all_data['all_species']['gene_changes'][gene] = {'snps': [0,0,0,0], 'all':[0,0,0,0], 'losses':[0,0,0,0], 'gains':[0,0,0,0]}
				# [obs, exp_between, exp_present, exp_pangeome]
			for change_type in change_types:	
				all_data['all_species']['gene_changes'][gene][change_type][0] += all_data[species_name]['gene_changes'][gene][change_type]

# Compute all species nulls by combining all trials across species:
null_types = ['between_host_genes', 'present_genes', 'pangenome_genes']
all_data['all_species']['null'] = {null_type : {} for null_type in null_types}

for species_name in all_data: 
	if species_name != 'all_species' : #CHANGE THIS
		for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
			species_null = all_data[species_name]['null'][null_type]
			all_species_null = all_data['all_species']['null'][null_type]
			for gene in species_null:
				if gene not in all_species_null:
					all_species_null[gene] = {'snps': [0]*num_trials, 'all': [0]*num_trials, 'gains': [0]*num_trials, 'losses': [0]*num_trials}
				for change_type in change_types:
					for i in range(0, num_trials):
						all_species_null[gene][change_type][i] += species_null[gene][change_type][i]

# Compute the expected number of changes under different nulls by averaging
all_gene_changes = all_data['all_species']['gene_changes']

for gene in all_gene_changes:
	for null_type in null_types:
		for change_type in change_types:
			expectation = numpy.array(all_data['all_species']['null'][null_type][gene][change_type]).mean() if gene in all_data['all_species']['null'][null_type] else 0
			if null_type == 'between_host_genes':
					all_gene_changes[gene][change_type][1] = expectation
			elif null_type == 'present_genes':
					all_gene_changes[gene][change_type][2] = expectation
			elif null_type == 'pangenome_genes':
					all_gene_changes[gene][change_type][3] = expectation

# Print the observed vs expected values of the genes in the all species gene changes:
gene_change_output = open('%s/gene_changes_across_species_%s.txt' %	(config.analysis_directory, cohort), 'w')

# Also store data for plotting cdf: (this says p-val, but actually expectations are being stored)
p_val_arrays={change_type: {null_type: [] for null_type in null_types} for change_type in change_types}

gene_order_sorted = sorted(all_gene_changes.items(), key=operator.itemgetter(1)) # Sorts by change_type indexed dictionary..?

for gene, _ in gene_order_sorted:
	string = gene
	for change_type in change_types:
		counts, p_val_between, p_val_present, p_val_pangenome = all_gene_changes[gene][change_type]
		# Fill in p_val_arrays with null expectations
		p_val_arrays[change_type]['between_host_genes'].append(p_val_between)		
		p_val_arrays[change_type]['present_genes'].append(p_val_present)		
		p_val_arrays[change_type]['pangenome_genes'].append(p_val_pangenome)
		# Print out observed and null number gene changes
		string += '\t' + str(counts) +'\t' +str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) 
	gene_change_output.write(string + '\n' )

gene_change_output.close()

##################################################
# list of common genes that show up							 #
##################################################

keywords={}
keywords['sugar'] = 'polysaccharide gluco sugar lact galac aldose glyco fructo manno xylan'.split()
keywords['fat'] = 'ester lip fatty enoyl-CoA'.split()
keywords['DNA'] = 'UspA DNA helicase deoxyribo chromosome replication excinuclease'.split()
keywords['RNA'] = 'RNA'.split()
keywords['transcription'] = 'transcript sigma ribosom'.split()
keywords['amino acid metabolism'] = 'glutamate dehydrogenase aminopeptidase chorismate'.split()
keywords['other metabolism'] = 'NADH ATP phosphohydrolase relaxase transferase hydrolase nucleotidase'.split()
keywords['reproduction'] = 'division'.split()
keywords['recombination'] = 'relaxosome recomb conjug transpos mobilization'.split() + ['mobile element']
keywords['antibiotic'] = 'HIPA'.split()
keywords['spore'] = 'sporulation spore'.split()
keywords['heat'] = 'heat'.split()
keywords['redox'] = 'aero redox reductase oxido FrrB'.split()
keywords['virus'] = 'phage integrase infection sialidase'.split()
keywords['membrane'] = 'membrane sortase chemotaxis HtpX flagella periplasm two-component wall capsul'.split()
keywords['transport'] = 'symporter transport channel secretion'.split() + ['sulfate permase']
keywords['metal'] = 'copper iron zinc metal'.split()
keywords['mucin'] = 'mucin'.split()
keywords['heparin'] = 'heparin'.split()
keywords['arginine'] = 'arginin'.split()
keywords['restriction-modification'] = 'restriction-modification'.split()
keywords['hypothetical'] = 'hypothetical unknown uncharacterized'.split()


# since this is a greedy algorithm, order the more important keywords first
keyword_order=['sugar', 'fat', 'DNA', 'RNA', 'transcription', 'amino acid metabolism', 'other metabolism', 'reproduction', 'recombination', 'antibiotic', 'spore', 'heat', 'redox', 'virus', 'membrane', 'transport', 'metal', 'mucin', 'heparin', 'arginine', 'restriction-modification', 'hypothetical']

common_genes={}

for keyword in keywords.keys():
		common_genes[keyword]={'snps': [0,0,0,0], 'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}

for gene in all_gene_changes:
	if all_gene_changes[gene]['all'][0] > 0 or all_gene_changes[gene]['snps'][0] > 0: 
		keyword_found=False
		for keyword in keyword_order:
			for regexp in keywords[keyword]:
				if regexp.lower() in gene.lower() and keyword_found==False:
					for change_type in change_types:
						for i in range(0,4):
							common_genes[keyword][change_type][i]+=all_gene_changes[gene][change_type][i]
					common_genes[keyword]['genes'].append(gene)
					keyword_found=True
		if keyword_found==False and gene !='':
			if gene not in common_genes.keys():
				common_genes[gene]={'snps': [0,0,0,0], 'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}
			for change_type in change_types:
				for i in range(0,4):
					common_genes[gene][change_type][i]+=all_gene_changes[gene][change_type][i]
			common_genes[gene]['genes'].append(gene)

#sorted list by total num

genes_sorted={}
for gene in common_genes.keys():
		genes_sorted[gene]=common_genes[gene]['all'][0]

sorted_genes = sorted(genes_sorted.items(), key=operator.itemgetter(1), reverse=True)

gene_change_keyword_output = open('%sgene_changes_across_species_keywords_%s.txt' % (config.analysis_directory, cohort),'w')
gene_change_keyword_output.write('keyword\tnum_all\texp_all_between\texp_all_present\texp_all_pangenome\tnum_gains\texp_gains_between\texp_gains_present\texp_gains_pangenome\tnum_loss\texp_loss_between\texp_loss_present\texp_loss_pangenome\tnum_snps\texp_snps_between\texp_snps_present\texp_snps_pangenome\tgene_names\n')

for i in range (0, len(sorted_genes)):
		gene=sorted_genes[i][0]		 
		string=gene
		for change_type in ['all', 'gains', 'losses', 'snps']: 
				for i in range(0,4):
						string += '\t' + str(common_genes[gene][change_type][i])
		string += '\t' + ';'.join(common_genes[gene]['genes'])
		gene_change_keyword_output.write(string+'\n')

gene_change_keyword_output.close()


###################################################
# plot CDF of p-values for the different nulls 
###################################################

pylab.figure(figsize=(6,6))
prevalence_axis = pylab.subplot(111)

prevalence_axis.set_ylabel('Fraction genes $\leq p$',labelpad=2)
prevalence_axis.set_xlabel('Expected number, $p$',labelpad=2)
prevalence_axis.set_xlim([0,20])
#prevalence_axis.set_ylim([0,1.1])

prevalence_axis.spines['top'].set_visible(False)
prevalence_axis.spines['right'].set_visible(False)
prevalence_axis.get_xaxis().tick_bottom()
prevalence_axis.get_yaxis().tick_left()

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.asarray(p_val_arrays['losses']['present_genes']))
prevalence_axis.step(xs,1-ns*1.0/ns[0],'b-',label='Within-host present genes',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.asarray(p_val_arrays['losses']['pangenome_genes']))
prevalence_axis.step(xs,1-ns*1.0/ns[0],'r-',label='Between-host gene changes',zorder=1)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.asarray(p_val_arrays['losses']['between_host_genes']))
prevalence_axis.step(xs,1-ns*1.0/ns[0],'k-',label='Pangenome',zorder=0)

#prevalence_axis.set_ylim([0,1.1])
#prevalence_axis.set_xlim([0,1.05])

prevalence_axis.legend(loc='upper right',frameon=False,fontsize=4)

pylab.savefig('%s/expected_gene_change_annotation_losses_%s.png' % (config.analysis_directory, cohort),bbox_inches='tight',dpi=300)
