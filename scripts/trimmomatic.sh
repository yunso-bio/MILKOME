#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N trim
#PBS -e trim.err
#PBS -o trim.log
#PBS -M yunso@dtu.dk
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=72:00:00

# Load modules
module load tools
module load perl
module load java/1.8.0
module load fastqc/0.11.9
module load trimmomatic/0.38

# Setup working directory and paths
path="path to working directory"
cd $path
file=$(cat "${path}/sample-list.txt")
raw="${path}/data/raw"
trimP="${path}/data/trimmomatic"
mkdir -p $trimP

# Run Trimommatic
for sample in $file; do
    f="$raw/${sample}_1.fq.gz"
    r="$raw/${sample}_2.fq.gz"
    echo 'processing' $sample
    java -jar /services/tools/trimmomatic/0.38/trimmomatic-0.38.jar PE -threads 40 -phred33 -trimlog logfile $f $r "${trimP}/${sample}_trimmed_1.fq.gz" "${trimP}/${sample}_trimmed_unpaired_1.fq.gz" "${trimP}/${sample}_trimmed_2.fq.gz" "${trimP}/${sample}_trimmed_unpaired_2.fq.gz" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done
