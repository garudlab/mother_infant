# ========================================================================
# Directory paths
# ========================================================================

data_dir=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/MegaMI/data
snps_dir=$data_dir/snps
genes_dir=$data_dir/genes
species_dir=$data_dir/species

scripts_dir=/u/project/ngarud/daisyche/mother_infant/scripts
shell_scripts_dir=$scripts_dir/shell_scripts

# ========================================================================
# Make snps/species_snps.txt: list of species which MIDAS called SNPS for
# + genes/species_genes.txt: species which MIDAS computed gene content for
# ========================================================================

ls $snps_dir > $snps_dir/species_snps.txt

# Manually remove non-species filenames from this list
# For this dataset, there are 502 species with SNP information

ls $genes_dir > $genes_dir/species_genes.txt

# Manually remove non-species filenames from this list
# For this dataset, there are 669 species with gene information
# Note that snps species are a strict subset of genes species

# ========================================================================
# Compress certain MIDAS output files to save space, code compatibility
# ========================================================================

# data/species
bzip2 $species_dir/count_reads.txt
bzip2 $species_dir/coverage.txt
bzip2 $species_dir/relative_abundance.txt
bzip2 $species_dir/species_prevalence.txt

# data/snps
qsub $shell_scripts_dir/qsub_compress_snps

# data/genes
qsub $shell_scripts_dir/qsub_compress_genes

# verify expected merged output after zipping
./$shell_scripts_dir/verify_midas_merge_output.sh

# ========================================================================
# Make a pickles folder
# ========================================================================

mkdir $data_dir/pickles

# ========================================================================
# Save species list text files under pickles (not really pickle but eh)
# ========================================================================

# Ran the following Python

from utils import parse_midas_data as pmd

pickles_dir = '/u/home/d/daisyche/dbd/data/pickles'

species_list = pmd.parse_depth_sorted_species_list()
file = open('%s/species_list.txt' % pickles_dir, 'w')
for spp in species_list:
    file.write(spp + '\n')

file.close()

good_species_list = pmd.parse_good_species_list()
file = open('%s/good_species_list.txt' % pickles_dir, 'w')

for spp in good_species_list:
    file.write(spp + '\n')

file.close()

# ========================================================================
# Postprocessing
# ========================================================================

python $scripts_dir/postprocess/calculate_core_genes.py

chmod +x $scripts_dir/postprocess/annotate_pvalue

qsub $scripts_dir/postprocess/qsub_postprocess_midas_data

$scripts_dir/shell_scripts/verify_midas_postprocess_output.sh

# ========================================================================
# Run all 'pickle' scripts, starting with pickle_qp.py,
# which should generate plotting script dependencies.
# ========================================================================

