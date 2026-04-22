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
mkdir -p ${path}/results/megahit
file=$(cat "path to a list of sampels")
trim_path="${path}/results/trimmomatic"

# Run MEGAHIT to assemble contigs
for sample in ${samples}; do
	f="${trim_path}/${sample}_trimmed_1.fq.gz"
	r="${trim_path}/${sample}_trimmed_2.fq.gz"
    echo 'processing' $file
    megahit -1 $f -2 $r -o "${path}/results/megahit/${sample}"
done

echo 'Job is done' | mail -s 'megahit is done' yunso@dtu.dk

exit
