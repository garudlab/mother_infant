# Constructs MegaMI MIDAS intermediate output database using symlinks

datadir=/u/project/ngarud/Garud_lab/metagenomic_fastq_files/

HMPdir=$datadir/HMP1-2
Backdir=$datadir/Backhed_2015
Ferrdir=$datadir/Ferretti_2018
Yassdir=$datadir/Yassour_2018
Shaodir=$datadir/Shao_2019
NIH1dir=$datadir/Olm_2019/NIH1
NIH2dir=$datadir/Olm_2019/NIH2
NIH3dir=$datadir/Olm_2019/NIH3
NIH4dir=$datadir/Olm_2019/NIH4
Sloan2dir=$datadir/Olm_2019/Sloan2

MegaMIdir=$datadir/MegaMI

# HMP1-2

samples=$(cat $HMPdir/HMP1-2_samples_final.txt)

for sample in $samples; do
	ln -s $HMPdir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Backhed_2015

samples=$(cat $Backdir/Backhed_run_accessions.txt)

for sample in $samples; do
	ln -s $Backdir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Ferretti_2018

samples=$(cat $Ferrdir/Ferretti_run_accessions.txt)

for sample in $samples; do
	ln -s $Ferrdir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Yassour_2018

samples=$(cat $Yassdir/Yassour_run_accessions.txt)

for sample in $samples; do
	ln -s $Yassdir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Shao_2019

samples=$(cat $Shaodir/PRJEB32631_run_accession_only.txt)

for sample in $samples; do
	ln -s $Shaodir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Olm_2019/NIH1

samples=$(cat $NIH1dir/Olm_NIH1_run_accessions_no_agg.txt)

for sample in $samples; do
	ln -s $NIH1dir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Olm_2019/NIH2

samples=$(cat $NIH2dir/Olm_NIH2_run_accessions_no_rep.txt)

for sample in $samples; do
	ln -s $NIH2dir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Olm_2019/NIH3

samples=$(cat $NIH3dir/Olm_NIH3_run_accessions_no_agg.txt)

for sample in $samples; do
	ln -s $NIH3dir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Olm_2019/NIH4

samples=$(cat $NIH4dir/Olm_NIH4_run_accessions_no_agg.txt)

for sample in $samples; do
	ln -s $NIH4dir/midas_output/$sample $MegaMIdir/midas_output/$sample
done

# Olm_2019/Sloan2

samples=$(cat $Sloan2dir/Olm_Sloan2_run_accessions_no_controls.txt)

for sample in $samples; do
	ln -s $Sloan2dir/midas_output/$sample $MegaMIdir/midas_output/$sample
done
