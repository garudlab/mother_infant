import sample_utils, config, parse_midas_data, sfs_utils, diversity_utils, gene_diversity_utils, core_gene_utils
import os, os.path, sys, gzip
import numpy

temporal_change_directory = '%s/temporal_changes/' % (config.data_directory)
intermediate_filename_template = '%s/%s.txt.gz'

min_coverage = config.min_median_coverage
min_sample_size = 2

def load_temporal_change_map(species_name, prev_cohort='all', min_coverage = 0):
		
		dir = "%s/cov%i_prev_%s" % (temporal_change_directory, min_coverage, prev_cohort)
		intermediate_filename = intermediate_filename_template % (dir, species_name)

		temporal_change_map = {}


		if not os.path.isfile(intermediate_filename):
				return temporal_change_map
		
		file = gzip.open(intermediate_filename,"r")
		file.readline() # header
		for line in file:
				items = line.split(",")
				if items[0].strip()!=species_name:
						continue
						
				sample_1 = items[1].strip()
				sample_2 = items[2].strip()
				type = items[3].strip()
				num_opportunities = float(items[4])
				perr = float(items[5])
				sample_pair = (sample_1, sample_2)
				if sample_pair not in temporal_change_map:
						temporal_change_map[sample_pair] = {}
				
				changes = []
				if len(items)<7:
						pass
				else:
						change_strs = items[6:]
						for change_str in change_strs:
						
								subitems = change_str.split(";")
								
								# switch on type of change
								if type=='snps':		
										gene_name = subitems[0].strip()
										contig = subitems[1].strip()
										position = long(subitems[2])
										variant_type = subitems[3].strip()
										A1 = float(subitems[4])
										D1 = float(subitems[5])
										A2 = float(subitems[6])
										D2 = float(subitems[7])
										changes.append( (gene_name, contig, position, variant_type, A1, D1, A2, D2) )
														
								elif type=='genes':
										gene_name = subitems[0].strip()
										D1 = float(subitems[1])
										Dm1 = float(subitems[2])
										D2 = float(subitems[3])
										Dm2 = float(subitems[4])
										changes.append( (gene_name, D1, Dm1, D2, Dm2) )
										
								elif type=='private_snps':
										
										gene_name = subitems[0].strip()
										contig = subitems[1].strip()
										position = long(subitems[2])
										variant_type = subitems[3].strip()
										A1 = float(subitems[4])
										D1 = float(subitems[5])
										A2 = float(subitems[6])
										D2 = float(subitems[7])
										changes.append( (gene_name, contig, position, variant_type, A1, D1, A2, D2) )
										
				temporal_change_map[sample_pair][type] = num_opportunities, perr, changes
		
		return temporal_change_map

def calculate_private_reversions_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):
		
		sample_pair = sample_1, sample_2
		if sample_pair not in temporal_change_map:
				return -1, None, None
				
		if 'private_snps' not in temporal_change_map[sample_pair]:
				return -1, None, None
				
		# otherwise, some hope! 
		
		private_snp_opportunities, private_snp_perr, private_snps = temporal_change_map[sample_pair]['private_snps']
		
		mutations = []
		private_snp_reversions = []
		for snp_change in private_snps:
		
				a,b,c,d,A1,D1,A2,D2 = snp_change
				
				if D1==0 or D2==0:
						private_snp_opportunities-=1
						continue
				
				f1 = A1*1.0/D1
				f2 = A2*1.0/D2
				
				if f1>=upper_threshold and f2<=lower_threshold:
						private_snp_reversions.append(snp_change)
				if f1<=upper_threshold and f2>=upper_threshold:
						mutations.append(snp_change)				
		
		return private_snp_opportunities, private_snp_perr, private_snp_reversions


def calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):

		sample_pair = sample_1, sample_2
		if sample_pair not in temporal_change_map:
				return -1, -1, [], []
				
		if 'snps' not in temporal_change_map[sample_pair]:
				return -1, -1, [], []
				
		# otherwise, some hope! 
		snp_opportunities, snp_perr, snp_changes = temporal_change_map[sample_pair]['snps']
		
		mutations = []
		reversions = []
		for snp_change in snp_changes:
		
				a,b,c,d,A1,D1,A2,D2 = snp_change
				
				f1 = A1*1.0/D1
				f2 = A2*1.0/D2
				
				if (f1<=lower_threshold) and (f2>=upper_threshold):
						mutations.append(snp_change)
				elif (f1>=upper_threshold) and (f2<=lower_threshold):
						reversions.append(snp_change)
						
		
		return snp_opportunities, snp_perr, mutations, reversions


def calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2, max_absent_copynum=config.gainloss_max_absent_copynum, min_normal_copynum=config.gainloss_min_normal_copynum, max_normal_copynum=config.gainloss_max_normal_copynum):


		sample_pair = sample_1, sample_2
		if sample_pair not in temporal_change_map:
				return -1, -1, [], []
				
		if 'genes' not in temporal_change_map[sample_pair]:
				return -1, -1, [], []
				
		# otherwise, some hope! 
		gene_opportunities, gene_perr, gene_changes = temporal_change_map[sample_pair]['genes']
		
		gains = []
		losses = []
		for gene_change in gene_changes:
		
				gene_name, D1, Dm1, D2, Dm2 = gene_change
				
				copynum_1 = D1/Dm1
				copynum_2 = D2/Dm2
				
				if (copynum_1<=max_absent_copynum) and (copynum_2>=min_normal_copynum) and (copynum_2<=max_normal_copynum):
						gains.append(gene_change)
				elif (copynum_2<=max_absent_copynum) and (copynum_1>=min_normal_copynum) and (copynum_1<=max_normal_copynum):
						losses.append(gene_change)
		
		return gene_opportunities, gene_perr, gains, losses
