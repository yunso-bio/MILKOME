#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N qc2
#PBS -e qc2.err
#PBS -o qc2.log
#PBS -M yunso@dtu.dk
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=12:00:00


module load tools
module load perl
module load java/1.8.0
module load fastqc/0.11.9
module load multiqc/1.12


cd /home/projects/dtu_00032/analysis/milk-cohort


file=$(cat "missing.txt")

mkdir qc/fastqc-host-removed
mkdir qc/multiqc-hr

for i in data/host_removed/*1.fastq.gz; do
    b=$(basename "$i" 1.fastq.gz)
    f="data/host_removed/${b}1.fastq.gz"
    r="data/host_removed/${b}2.fastq.gz"
	fastqc -o qc/fastqc-host-removed -t 40 $f $r
done

echo
echo
echo

n=$(ls qc/fastqc-host-removed/*1_fastqc.html | wc -l)
echo $n

multiqc qc/fastqc-host-removed multiqc-HR


exit
