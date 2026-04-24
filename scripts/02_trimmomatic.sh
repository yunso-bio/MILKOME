#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N trim
#PBS -e trim.err
#PBS -o trim.log
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
BASE="path to working directory"
samples=$(cat "Path to a list of samples")
raw_path="$BASE/data/raw"
trim_path="$BASE/results/trimmomatic"
mkdir -p $trim_path

# Run Trimommatic
for sample in ${samples}; do
    f="$raw_path/${sample}_1.fq.gz"
    r="$raw_path/${sample}_2.fq.gz"
    echo 'processing' $sample
    java -jar /services/tools/trimmomatic/0.38/trimmomatic-0.38.jar PE -threads 40 -phred33 "$trim_path/${sample}_trimlog.txt" $f $r "$trim_path/${sample}_trimmed_1.fq.gz" "$trim_path/${sample}_trimmed_unpaired_1.fq.gz" "$trim_path/${sample}_trimmed_2.fq.gz" "$trim_path/${sample}_trimmed_unpaired_2.fq.gz" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
done
