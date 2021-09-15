#!/bin/bash

# ===================================================================================
# Verifies that MIDAS postprocess output is complete
# 
# Specifically, this script checks that:
# 	- for each species in snps, there is (see snps_files below)
# 	- for each species in genes, there is (see genes_files below)
# ===================================================================================

proj_dir=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/MegaMI

species_snps=$(cat $proj_dir/data/snps/species_snps.txt)
species_genes=$(cat $proj_dir/data/genes/species_genes.txt)

snps_dir=$proj_dir/data/snps
genes_dir=$proj_dir/data/genes

snps_files='readme.txt  snps_alt_allele.txt.bz2  snps_depth.txt.bz2  snps_info.txt.bz2  snps_log.txt  snps_ref_freq.txt.bz2  snps_summary.txt annotated_snps.txt.bz2 coverage_distribution.txt.bz2 full_coverage_distribution.txt.bz2 gene_coverage.txt.bz2 marker_coverage.txt.bz2 within_sample_sfs.txt.bz2'
unzipped_snps_files='snps_alt_allele.txt  snps_depth.txt  snps_info.txt  snps_ref_freq.txt'

genes_files='genes_copynum.txt.bz2  genes_depth.txt.bz2  genes_presabs.txt.bz2  genes_reads.txt.bz2  genes_summary.txt  readme.txt'
unzipped_genes_files='genes_copynum.txt  genes_depth.txt  genes_presabs.txt  genes_reads.txt'

# Check snps

for species in $species_snps; do
	
	species_dir=$snps_dir/$species
	
	for file in $snps_files; do
		if [ ! -f $species_dir/$file ];
		then echo $species is missing snps information; break; fi
	done
	
	for file in $unzipped_snps_files; do
		if [ -f $species_dir/$file ];
		then echo $species has unzipped snps files; break; fi
	done
	
done

# Check genes

for species in $species_genes; do
	
	species_dir=$genes_dir/$species
	
	for file in $genes_files; do
		if [ ! -f $species_dir/$file ];
		then echo $species is missing genes information; break; fi
	done
	
	for file in $unzipped_genes_files; do
		if [ -f $species_dir/$file ];
		then echo $species has unzipped genes files; break; fi
	done
	
done

