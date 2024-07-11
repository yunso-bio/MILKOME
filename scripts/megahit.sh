#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N megahit
#PBS -e megahit.err
#PBS -o megahit.log
#PBS -M yunso@dtu.dk
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=72:00:00

cd /home/projects/dtu_00032/analysis/milk-cohort

module load tools
module load megahit/1.2.9

mkdir megahit

file=$(cat "10-firsthalf.txt")
for i in $file; do
	f=$"/home/projects/dtu_00032/analysis/milk-cohort/data/trimmomatic/${i}_trimmed_1.fq.gz"
	r=$"/home/projects/dtu_00032/analysis/milk-cohort/data/trimmomatic/${i}_trimmed_2.fq.gz"
    echo 'processing' $file
    megahit -1 $f -2 $r -o "data/megahit/${i}"
done


echo 'Job is done' | mail -s 'megahit is done' yunso@dtu.dk

exit
