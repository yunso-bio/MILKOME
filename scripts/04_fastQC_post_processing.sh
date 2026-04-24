#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N qc2
#PBS -e qc2.err
#PBS -o qc2.log
#PBS -M yunso@dtu.dk
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=12:00:00

# Load modules
module load tools
module load perl
module load java/1.8.0
module load fastqc/0.11.9
module load multiqc/1.12

# Setup working directory and paths
BASE="path to working directory"
mkdir -p ${path}/results/qc/{fastqc_cleaned,multiqc_cleaned}
stool_samples=$(cat "Path to stool sample list")
enriched_samples=$(cat "Path to enriched sample list")

# Run FastQC for stool samples
for sample in ${stool_samples}; do
    f="$BASE/results/host_removed/${sample}_1.fq.gz"
    r="$BASE/results/host_removed/${sample}_2.fq.gz"
	fastqc -o ${path}/results/qc/fastqc_cleaned -t 40 $f $r
done

# Run FastQC for enriched samples
for sample in ${enriched_samples}; do
    f="$BASE/results/trimmomatic/${sample}_trimmed_1.fq.gz"
    r="$BASE/results/trimmomatic/${sample}_trimmed_2.fq.gz"
	fastqc -o ${path}/results/qc/fastqc -t 40 $f $r
done

n=$(ls qc/fastqc_cleaned/*1_fastqc.html | wc -l)
echo $n

# Run MultiQC
multiqc qc/fastqc_cleaned multiqc_cleaned

exit
