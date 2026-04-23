#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N qc
#PBS -e qc.err
#PBS -o qc.log
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
path="path to working directory"
mkdir -p ${path}/results/qc/{fastqc,multiqc}
samples=$(cat "Path to sample list file")

# Run FastQC
for sample in ${samples}; do
    f="$path/data/raw/${sample}_1.fq.gz"
    r="$path/data/raw/${sample}_2.fq.gz"
	fastqc -o ${path}/results/qc/fastqc -t 40 $f $r
done

# print out the number of processed samples
echo
n=$(ls $path/results/qc/fastqc/*1_fastqc.html | wc -l)
echo $n

# Run multiQC
multiqc -o "$path/results/qc/multiqc" "$path/results/qc/fastqc"

exit
