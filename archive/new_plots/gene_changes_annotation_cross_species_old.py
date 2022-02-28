# Within-host gene change annotation

import os
import sys
import numpy
import random
import pickle
import operator

from utils import sample_utils as su, config, parse_midas_data

#######################################
# code for loading cross-species data #
#######################################
num_trials=100

good_species_list = parse_midas_data.parse_good_species_list() 
cohorts = ['backhed', 'ferretti', 'yassour', 'shao', 'olm', 'hmp']

ddir = config.data_directory
pdir = "%s/pickles" % ddir

all_data={}
for species_name in good_species_list:
		print species_name
		if os.path.isfile("%s/gene_changes_anno_old_%s.pkl" % (pdir, species_name)):
				all_data_species = pickle.load( open( "%s/gene_changes_anno_old_%s.pkl" % (pdir, species_name), "rb" ))
				if (len(all_data_species['hmp'][species_name].keys()) >0): # check if there were any gene changes to be outputted. 
						new_all_data_speces = {}
						for cohort in all_data_species.keys():
							new_all_data_speces[cohort] = all_data_species[cohort][species_name]
						all_data[species_name] = new_all_data_speces

#	 sum	all gene changes	across species
all_data['all_species']={cohort: {} for cohort in cohorts}
#for gene_type in ['gene_changes','gene_changes_category']:
for gene_type in ['gene_changes']:
		for cohort in cohorts:
			all_data['all_species'][cohort][gene_type]={}
			#all_data['all_species']['gene_changes_category']={}

for species_name in all_data.keys():
		print species_name
		if species_name != 'all_species':
			for cohort in all_data[species_name]:
				for gene_type in ['gene_changes']:
						if gene_type in all_data[species_name][cohort]:
							for gene in all_data[species_name][cohort][gene_type].keys():
									if gene not in all_data['all_species'][cohort][gene_type].keys():
											all_data['all_species'][cohort][gene_type][gene] = {'all':[0,0,0,0],'losses':[0,0,0,0],'gains':[0,0,0,0]} #[obs, exp_between, exp_present, exp_pangeome]
									for change_type in ['gains','losses','all']:	
											all_data['all_species'][cohort][gene_type][gene][change_type][0] += all_data[species_name][cohort][gene_type][gene][change_type]
																																 

# compute a null for all species by combining all trials across species:
for cohort in cohorts:
	all_data['all_species'][cohort]['null']={}
	#for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes', 'between_host_category','present_category','pangenome_category']:
	for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
			all_data['all_species'][cohort]['null'][null_type]={}

for species_name in all_data.keys(): 
		print species_name
		for cohort in all_data[species_name]:
			if 'null' in all_data[species_name][cohort]:
				if species_name != 'all_species' : #CHANGE THIS
						for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
								for gene in all_data[species_name][cohort]['null'][null_type].keys():
										if gene not in	all_data['all_species'][cohort]['null'][null_type].keys():
												all_data['all_species'][cohort]['null'][null_type][gene]={'all':[0]*num_trials,'gains':[0]*num_trials,'losses':[0]*num_trials}
										for change_type in ['gains','losses','all']:
												for i in range(0, num_trials):
														all_data['all_species'][cohort]['null'][null_type][gene][change_type][i]+=all_data[species_name][cohort]['null'][null_type][gene][change_type][i]

# compute the expected number of changes under different nulls by averaging
for cohort in cohorts:
	for gene in all_data['all_species'][cohort]['gene_changes'].keys():
			for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
					for change_type in ['gains','losses','all']:
							if gene in all_data['all_species'][cohort]['null'][null_type].keys():
									expectation = numpy.array(all_data['all_species'][cohort]['null'][null_type][gene][change_type]).mean()
									#prob=sum(all_data['all_species']['gene_changes'][gene][change_type] <= tmp_null)/float(num_trials)
							else:
									expectation=0
							if null_type=='between_host_genes':
									all_data['all_species'][cohort]['gene_changes'][gene][change_type][1]=expectation
							elif null_type=='present_genes':
									all_data['all_species'][cohort]['gene_changes'][gene][change_type][2]=expectation
							else:
									all_data['all_species'][cohort]['gene_changes'][gene][change_type][3]=expectation		

'''
# repeat for gene_category

for gene in all_data['all_species']['gene_changes_category'].keys():
		for null_type in ['between_host_category','present_category','pangenome_category']:
				if gene in all_data['all_species']['null'][null_type].keys():
						tmp_null = numpy.array(all_data['all_species']['null'][null_type][gene])
						prob=sum(all_data['all_species']['gene_changes_category'][gene][0] <= tmp_null)/float(num_trials)
				else:
						prob=0
				if null_type=='between_host_category':
						all_data['all_species']['gene_changes_category'][gene][1]=prob
				elif null_type=='present_category':
						all_data['all_species']['gene_changes_category'][gene][2]=prob
				else:
						all_data['all_species']['gene_changes_category'][gene][3]=prob
				
'''

for cohort in cohorts:
	
	# print the observed vs expected values of the genes in the all species gene changes:
	outFile_gene_change=open('%sgene_changes_across_species_v2_%s.txt' %	 (config.analysis_directory, cohort) ,'w')
	
	# store data for plotting cdf: (this says p-val, but actually expectations are being stored)
	p_val_arrays={}
	for change_type in ['all','gains','losses']: 
			p_val_arrays[change_type]={'between':[],'present':[],'pangenome':[]}
	
	gene_order_sorted = sorted(all_data['all_species'][cohort]['gene_changes'].items(), key=operator.itemgetter(1))
	for i in range(0, len(gene_order_sorted)):
			gene=gene_order_sorted[i][0]
			string=gene 
			for change_type in ['all','gains','losses']:
					counts=all_data['all_species'][cohort]['gene_changes'][gene][change_type][0]
					p_val_between=all_data['all_species'][cohort]['gene_changes'][gene][change_type][1]
					p_val_arrays[change_type]['between'].append(p_val_between)
					p_val_present=all_data['all_species'][cohort]['gene_changes'][gene][change_type][2]
					p_val_arrays[change_type]['present'].append(p_val_present)
					p_val_pangenome=all_data['all_species'][cohort]['gene_changes'][gene][change_type][3]
					p_val_arrays[change_type]['pangenome'].append(p_val_pangenome)
					string += '\t' + str(counts) +'\t' +str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) 
			print string
			outFile_gene_change.write(string + '\n' )
	
	outFile_gene_change.close()



##################################################
# list of common genes that show up							 #
##################################################

keywords={}
keywords['ABC transporter']=['ABC']
keywords['phage']=['hage']
keywords['transposon']=['onjugati','anspos']
keywords['mobilization']=['mob','obilization','obile']
keywords['integrase']=['ntegrase']
keywords['plasmid']=['plasmid']
keywords['recombinase']=['ecombinase']
keywords['tRNA']=['tRNA']
keywords['ATP']=['ATP']
keywords['excisionase']=['xcisionase']
keywords['transmembrane']=['embrane']
keywords['replication']=['eplication']
keywords['regulator']=['egulator']
keywords['transcription']=['anscription']
keywords['toxin']=['toxin']
keywords['restriction']=['estriction']
keywords['replication']=['eplication']
keywords['transferase']=['ansferase']
keywords['reductase']=['eductase']
keywords['phosphatase']=['phosphatase']
keywords['helicase']=['elicase']
keywords['kinase']=['kinase']
keywords['dehydrogenase']=['dehydrogenase']
keywords['drug']=['drug']
keywords['cell wall']=['ell wall']
keywords['primase']=['imase']
keywords['resistance']=['resistance']
keywords['hydrolase']=['ydrolase']
keywords['topoisomerase']=['opoisomerase']
keywords['hypothetical'] = ['ypothetical']
#keywords['other']=['other']

# since this is a greedy algorithm, order the more important keywords first
keyword_order=['ABC transporter','phage','transposon','mobilization','integrase', 'plasmid','recombinase','tRNA','ATP','excisionase','transmembrane','replication','regulator','transcription','toxin','restriction','replication','transferase','reductase','phosphatase','helicase','kinase','dehydrogenase','drug','cell wall','primase','resistance','hydrolase','topoisomerase','hypothetical']

for cohort in cohorts:
	common_genes={}
	
	# key=keyword
	# value={}
	# num={}
	# genes={}
	for keyword in keywords.keys():
			common_genes[keyword]={'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}
	
	for gene in all_data['all_species'][cohort]['gene_changes']:
		if all_data['all_species'][cohort]['gene_changes'][gene]['all'][0] > 0: 
			keyword_found=False
			for keyword in keyword_order:
					for regexp in keywords[keyword]:
							if regexp in gene and keyword_found==False:
									for change_type in ['all','gains','losses']:
											for i in range(0,4):
													common_genes[keyword][change_type][i]+=all_data['all_species'][cohort]['gene_changes'][gene][change_type][i]
									common_genes[keyword]['genes'].append(gene)
									keyword_found=True
			if keyword_found==False and gene !='':
					if gene not in common_genes.keys():
							common_genes[gene]={'all':[0,0,0,0], 'gains':[0,0,0,0], 'losses':[0,0,0,0], 'genes':[]}
					for change_type in ['all','gains','losses']:
							for i in range(0,4):
									common_genes[gene][change_type][i]+=all_data['all_species'][cohort]['gene_changes'][gene][change_type][i]
					common_genes[gene]['genes'].append(gene)
	
	genes_sorted={}
	for gene in common_genes.keys():
			genes_sorted[gene]=common_genes[gene]['all'][0]
	
	sorted_genes = sorted(genes_sorted.items(), key=operator.itemgetter(1), reverse=True)
	
	outFile_keywords=open('%s/gene_changes_across_species_keywords_v2_%s.txt' %	 (config.analysis_directory, cohort),'w')
	
	outFile_keywords.write('keyword\tnum_all\texp_all_between\texp_all_present\texp_all_pangenome\tnum_gains\texp_gains_between\texp_gains_present\texp_gains_pangenome\tnum_loss\texp_loss_between\texp_loss_present\texp_loss_pangenome\tgene_names\n')
	
	for i in range (0, len(sorted_genes)):
			gene=sorted_genes[i][0]		 
			string=gene
			for change_type in ['all','gains','losses']: 
					for i in range(0,4):
							string += '\t' + str(common_genes[gene][change_type][i])
			string += '\t' + ';'.join(common_genes[gene]['genes'])
			print string
			outFile_keywords.write(string+'\n')
	
	outFile_keywords.close()
