#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N megahit
#PBS -e megahit.err
#PBS -o megahit.log
#PBS -M yunso@dtu.dk
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=72:00:00

# Load modules 
module load tools
module load megahit/1.2.9

# Setup working directory and paths
path="path to working directory"
mkdir -p "$path/results/megahit"

stool_samples=$(cat "path to a list of stool samples")
clean_path="$path/results/host_removed"

# Run MEGAHIT for stool samples to assemble contigs
for sample in ${stool_samples}; do
    f="$clean_path/${sample}_1.fq.gz"
    r="$clean_path/${sample}_2.fq.gz"
    echo "processing $sample"
    megahit -1 "$f" -2 "$r" -o "$path/results/megahit/$sample"
done

# Setup working directory and paths
enriched_samples=$(cat "path to a list of enriched samples")
trim_path="$path/results/trimmomatic"

# Run MEGAHIT for enriched samples to assemble contigs
for sample in ${enriched_samples}; do
    f="$trim_path/${sample}_trimmed_1.fq.gz"
    r="$trim_path/${sample}_trimmed_2.fq.gz"
    echo "processing $sample"
    megahit -1 "$f" -2 "$r" -o "$path/results/megahit/$sample"
done

echo 'Job is done' | mail -s 'megahit is done' yunso@dtu.dk

exit
