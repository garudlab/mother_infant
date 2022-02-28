# This uses the postprocessed data files
import matplotlib
matplotlib.use('Agg') 
import parse_midas_data
import sample_utils
import sys
import numpy
import diversity_utils
import stats_utils
from math import log10,ceil
import matplotlib.pyplot as plt
import config
import sfs_utils
import os.path


minfreq=1
upper_threshold = minfreq
lower_threshold = 1 - upper_threshold
min_coverage = 20

good_species_list=parse_midas_data.parse_good_species_list()


def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color, linestyle='-')
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)    
    plt.setp(bp['fliers'], color=color, markersize=3)


# set up dataframe to store data
dataset_cohort={}
dataset_cohort['backhed']=['B','4M','12M','M']
dataset_cohort['yassour']=['MGest', 'MBirth', 'M3', 'CBirth', 'C14', 'C1', 'C2', 'C3']
dataset_cohort['ferretti']=['M0','I1','I2','I3','I4','I5']

ref_freq_frac_QP={}
ref_freq_frac_nonQP={}
within_sites_1D_species_QP={}
total_sites_1D_species_QP={}
within_sites_4D_species_QP={}
total_sites_4D_species_QP={}
within_sites_1D_species_nonQP={}
total_sites_1D_species_nonQP={}
within_sites_4D_species_nonQP={}
total_sites_4D_species_nonQP={}
for dataset in dataset_cohort:     
    for age in dataset_cohort[dataset]:   
        ref_freq_frac_QP[age]=numpy.array([])  
        ref_freq_frac_nonQP[age]=numpy.array([])
        within_sites_1D_species_QP[age]=0
        total_sites_1D_species_QP[age]=0
        within_sites_4D_species_QP[age]=0
        total_sites_4D_species_QP[age]=0
        within_sites_1D_species_nonQP[age]=0
        total_sites_1D_species_nonQP[age]=0
        within_sites_4D_species_nonQP[age]=0
        total_sites_4D_species_nonQP[age]=0

for species_name in good_species_list:
  #
  # Load subject and sample metadata
  sys.stderr.write("Loading sample metadata...\n")
  subject_sample_map = sample_utils.parse_subject_sample_map()
  sample_order_map = sample_utils.parse_sample_order_map()
  sys.stderr.write("Done!\n") 
  #
  # Load SNP information for species_name
  sys.stderr.write("Loading SFSs for %s...\t" % species_name)
  samples, sfs_map_4D = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
  samples_1D, sfs_map_1D = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D'])) 
  sys.stderr.write("Done!\n")
  #
  # Load genomic coverage distributions
  sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
  median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
  sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
  samples = numpy.array(samples)
  samples_QP = diversity_utils.calculate_haploid_samples(species_name)
  # make array of non-QP samples:
  samples_nonQP=numpy.array([])
  for sample in samples:
      if sample not in samples_QP:
          samples_nonQP=numpy.append(samples_nonQP, sample)
  #
  median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])
  #
  # Only plot samples above a certain depth threshold
  desired_samples = samples[(median_coverages>=min_coverage)]
  desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])

  
  if len(desired_samples) > 0:
      #
      ###################################
      #
      # Calculate within polymorphism rates
      #
      ###################################
      #
      sample_names = []
      #
      between_rates_1D = []
      between_rates_4D = []
      #
      within_rates_1D = []
      within_rates_4D = []
      within_rate_lowers = []
      within_rate_uppers = []
      #
      median_depths = []
      depth_lowers = []
      depth_uppers = []
      #
      within_sites_1D_array=[]
      total_sites_1D_array=[]
      within_sites_4D_array=[]
      total_sites_4D_array=[]
      #
      for sample in desired_samples:
          # plot pi (or, num intermediate alleles)
          #within_sites_1D, between_sites_1D, total_sites_1D = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map_1D[sample], lower_threshold, upper_threshold)
          #within_sites_4D, between_sites_4D, total_sites_4D = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map_4D[sample], lower_threshold, upper_threshold)
          # plot s (or, num low freq segregating sites)
          within_sites_1D, between_sites_1D, total_sites_1D = sfs_utils.calculate_singleton_rates_from_sfs_map(sfs_map_1D[sample], lower_threshold=0, upper_threshold=0.05)
          within_sites_4D, between_sites_4D, total_sites_4D = sfs_utils.calculate_singleton_rates_from_sfs_map(sfs_map_4D[sample], lower_threshold=0, upper_threshold=0.05)
          if total_sites_1D >0 and total_sites_4D >0:  
              within_rate_1D = within_sites_1D*1.0/total_sites_1D
              between_rate_1D = between_sites_1D*1.0/total_sites_1D
              between_rates_1D.append(between_rate_1D)
              within_rates_1D.append(within_rate_1D)
              #
              within_sites_1D_array.append(within_sites_1D)
              total_sites_1D_array.append(total_sites_1D)
              #
              within_rate_4D = within_sites_4D*1.0/total_sites_4D
              between_rate_4D = between_sites_4D*1.0/total_sites_4D
              between_rates_4D.append(between_rate_4D)
              within_rates_4D.append(within_rate_4D)
              #
              within_sites_4D_array.append(within_sites_4D)
              total_sites_4D_array.append(total_sites_4D)              
              #
              #
              # Calculate 95% confidence intervals
              #within_rate_lower_4D, within_rate_upper_4D = stats_utils.calculate_poisson_rate_interval(within_sites_4D, total_sites_1D,alpha=0.05)
              #within_rate_lowers_4D.append(within_rate_lower_4D)
              #within_rate_uppers_4D.append(within_rate_upper_4D)
              # 
              #depths, counts = sfs_utils.calculate_depth_distribution_from_sfs_map(sfs_map_4D[sample])
              #dlower, dupper = stats_utils.calculate_IQR_from_distribution(depths, counts)
              #dmedian = stats_utils.calculate_median_from_distribution(depths,counts)
              #
              #depth_lowers.append(dlower)
              #depth_uppers.append(dupper)
              #median_depths.append(dmedian)
              sample_names.append(sample)

      # Sort them all in descending order of within-host diversity		
      within_rates_4D, between_rates_4D, within_sites_4D_array, total_sites_4D_array, within_rates_1D, between_rates_1D, within_sites_1D_array, total_sites_1D_array,  sample_names = (numpy.array(x) for x in zip(*sorted(zip(within_rates_4D, between_rates_4D, within_sites_4D_array, total_sites_4D_array, within_rates_1D, between_rates_1D, within_sites_1D_array, total_sites_1D_array,  sample_names),reverse=True)))

      #within_rate_lowers = numpy.clip(within_rate_lowers, 1e-09,1)



  ###################################
  #
  # Partition data by age range
  #
  ###################################

  # order that we want:

  #CBirth Yassour Birth -- what does that mean? 
  #I1 Ferretti 1 day

  #I2 Ferretti 3 days

  #B Backhed 4 days

  #I3 Ferretti 7 days

  #C14 Yassour 14 days

  #C1 Yassour 1 month
  #I4 Ferretti 1 month

  #C2 Yassour 2 months

  #C3 Yassour 3 months

  #4M Backhed 4 months
  #I5 Ferretti 4 months

  #12M Backhed 1 year

  #M Backhed Mother
  #MGest Yassour
  #MBirth Yassour
  #M3 Yassour
  #M0 Ferretti

  #HMP?

  #store refreq for each age for QP samples:
  for dataset in dataset_cohort:
          for age in dataset_cohort[dataset]:
                  sample_names_cohort=[]
                  sample_names_cohort = sample_utils.get_sample_names(dataset, age)
                  for sample_name in sample_names_cohort:
                      if sample_name in sample_names:
                          idx=numpy.where(sample_names==sample_name)
                          if sample_name in samples_QP:
                              ref_freq_frac_QP[age]= numpy.append(ref_freq_frac_QP[age], within_rates_4D[idx])
                              within_sites_1D_species_QP[age]+=within_sites_1D_array[idx]
                              total_sites_1D_species_QP[age]+=total_sites_1D_array[idx]
                              within_sites_4D_species_QP[age]+=within_sites_4D_array[idx]
                              total_sites_4D_species_QP[age]+=total_sites_4D_array[idx]
                          elif sample_name in samples_nonQP:
                              ref_freq_frac_nonQP[age]= numpy.append(ref_freq_frac_nonQP[age], within_rates_4D[idx])
                              within_sites_1D_species_nonQP[age]+=within_sites_1D_array[idx]
                              total_sites_1D_species_nonQP[age]+=total_sites_1D_array[idx]
                              within_sites_4D_species_nonQP[age]+=within_sites_4D_array[idx]
                              total_sites_4D_species_nonQP[age]+=total_sites_4D_array[idx]



######################################################
# after iterating through all species, now plot
######################################################

# combine all 1 month samples
ref_freq_frac_QP['1_month']=numpy.concatenate((ref_freq_frac_QP['C1'],ref_freq_frac_QP['I4']), axis=0)
ref_freq_frac_nonQP['1_month']=numpy.concatenate((ref_freq_frac_nonQP['C1'],ref_freq_frac_nonQP['I4']), axis=0)

within_sites_1D_species_QP['1_month'] = within_sites_1D_species_QP['C1'] + within_sites_1D_species_QP['I4']
total_sites_1D_species_QP['1_month'] = total_sites_1D_species_QP['C1'] + total_sites_1D_species_QP['I4']
within_sites_4D_species_QP['1_month'] = within_sites_4D_species_QP['C1'] + within_sites_4D_species_QP['I4']
total_sites_4D_species_QP['1_month'] = total_sites_4D_species_QP['C1'] + total_sites_4D_species_QP['I4']

within_sites_1D_species_nonQP['1_month'] = within_sites_1D_species_nonQP['C1'] + within_sites_1D_species_nonQP['I4']
total_sites_1D_species_nonQP['1_month'] = total_sites_1D_species_nonQP['C1'] + total_sites_1D_species_nonQP['I4']
within_sites_4D_species_nonQP['1_month'] = within_sites_4D_species_nonQP['C1'] + within_sites_4D_species_nonQP['I4']
total_sites_4D_species_nonQP['1_month'] = total_sites_4D_species_nonQP['C1'] + total_sites_4D_species_nonQP['I4']


# combine all 4 month sample
ref_freq_frac_QP['4_months']=numpy.concatenate((ref_freq_frac_QP['4M'],ref_freq_frac_QP['I5']), axis=0)
ref_freq_frac_nonQP['4_months']=numpy.concatenate((ref_freq_frac_nonQP['4M'],ref_freq_frac_nonQP['I5']), axis=0)


within_sites_1D_species_QP['4_months'] = within_sites_1D_species_QP['4M'] + within_sites_1D_species_QP['I5']
total_sites_1D_species_QP['4_months'] = total_sites_1D_species_QP['4M'] + total_sites_1D_species_QP['I5']
within_sites_4D_species_QP['4_months'] = within_sites_4D_species_QP['4M'] + within_sites_4D_species_QP['I5']
total_sites_4D_species_QP['4_months'] = total_sites_4D_species_QP['4M'] + total_sites_4D_species_QP['I5']

within_sites_1D_species_nonQP['4_months'] = within_sites_1D_species_nonQP['4M'] + within_sites_1D_species_nonQP['I5']
total_sites_1D_species_nonQP['4_months'] = total_sites_1D_species_nonQP['4M'] + total_sites_1D_species_nonQP['I5']
within_sites_4D_species_nonQP['4_months'] = within_sites_4D_species_nonQP['4M']+ within_sites_4D_species_nonQP['I5']
total_sites_4D_species_nonQP['4_months'] = total_sites_4D_species_nonQP['4M']+total_sites_4D_species_nonQP['I5']


# combine all mother samples
ref_freq_frac_QP['all_mothers']=numpy.concatenate((ref_freq_frac_QP['M'],ref_freq_frac_QP['MBirth'],ref_freq_frac_QP['M0']), axis=0)
ref_freq_frac_nonQP['all_mothers']=numpy.concatenate((ref_freq_frac_nonQP['M'],ref_freq_frac_nonQP['MBirth'],ref_freq_frac_nonQP['M0']), axis=0)


within_sites_1D_species_QP['all_mothers'] = within_sites_1D_species_QP['M']+ within_sites_1D_species_QP['MBirth'] + within_sites_1D_species_QP['M0']
total_sites_1D_species_QP['all_mothers'] = total_sites_1D_species_QP['M']+ total_sites_1D_species_QP['MBirth']+ total_sites_1D_species_QP['M0']
within_sites_4D_species_QP['all_mothers'] = within_sites_4D_species_QP['M']+within_sites_4D_species_QP['MBirth']+ within_sites_4D_species_QP['M0']
total_sites_4D_species_QP['all_mothers'] = total_sites_4D_species_QP['M']+ total_sites_4D_species_QP['MBirth']+ total_sites_4D_species_QP['M0']

within_sites_1D_species_nonQP['all_mothers'] = within_sites_1D_species_nonQP['M']+within_sites_1D_species_nonQP['MBirth']+ within_sites_1D_species_nonQP['M0']
total_sites_1D_species_nonQP['all_mothers'] = total_sites_1D_species_nonQP['M']+total_sites_1D_species_nonQP['MBirth']+ total_sites_1D_species_nonQP['M0']
within_sites_4D_species_nonQP['all_mothers'] = within_sites_4D_species_nonQP['M']+within_sites_4D_species_nonQP['MBirth']+ total_sites_1D_species_nonQP['M0']
total_sites_4D_species_nonQP['all_mothers'] = total_sites_4D_species_nonQP['M']+total_sites_4D_species_nonQP['MBirth']+ total_sites_1D_species_nonQP['M0']



age_order=['CBirth','I1','I2','B','I3', 'C14', '1_month','C2','C3','4_months','12M','all_mothers']


colors=['#fed976','#feb24c', '#fd8d3c','#fc4e2a', '#e31a1c', '#b10026', '#bdd7e7','#6baed6', '#08519c','#8c510a','pink','#d9ef8b', '#66bd63', 'black','red','green','blue','yellow']

####################################
# Compute pn/ps for each age class #
####################################


outfile=open('pnps.txt', 'w')
pnps={}
pnps['QP'] = {}
pnps['nonQP'] = {}   
pss=[]
pnps_QP=[]
pnps_nonQP=[]
for age in  age_order: 
    pn=within_sites_1D_species_nonQP[age]*1.0/(total_sites_1D_species_nonQP[age] +1.0)
    ps=(within_sites_4D_species_nonQP[age]*1.0+1.0)/(total_sites_4D_species_nonQP[age] + 1.0)
    pnps['nonQP'][age] = pn/ps
    pn=(within_sites_1D_species_QP[age]*1.0+1.0)/(total_sites_1D_species_QP[age] +1.0)
    ps=(within_sites_4D_species_QP[age]*1.0+1.0)/(total_sites_4D_species_QP[age] + 1.0)
    pnps['QP'][age] = pn/ps
    outfile.write(age + '\t' + str(ps) + '\t' + str(pnps['QP'][age]) + '\t' + str(pnps['nonQP'][age]) + '\n' )
    pss.append(ps)
    pnps_QP.append(pnps['QP'][age])
    pnps_nonQP.append(pnps['nonQP'][age])  
outfile.close()


# plot

fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 6
fig_size[1] = 6

plt.plot(pss, pnps_QP, 'go--', linewidth=2, markersize=12)
plt.plot(pss, pnps_nonQP, 'bo--', linewidth=2, markersize=12)

plt.savefig(os.path.expanduser('~/mother_infant/analysis/pnps_ps.png'), dpi=600)
