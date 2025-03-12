#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N mouse-gene-removal
#PBS -e mouse-gene-removal.err
#PBS -o mouse-gene-removal.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=48:00:00


module load tools
module load bowtie2/2.5.3

cd /home/projects/dtu_00032/analysis/milk-cohort/mouse

rm -r data/host_removed
rm -r data/bowtie
mkdir data/host_removed
mkdir data/bowtie

file=$(cat "/home/projects/dtu_00032/analysis/milk-cohort/sample-list.txt")

for N in $file; do
        f=$"/home/projects/dtu_00032/analysis/milk-cohort/mouse/data/host_removed/${N}.1"
        r=$"/home/projects/dtu_00032/analysis/milk-cohort/mouse/data/host_removed/${N}.2"
        echo 'processing' ${N}
        echo 
        echo 'processing' $N
        bowtie2 -p 40 -x /home/projects/dtu_00032/db/GRCm39/GRCm39 -1 $f -2 $r --un-conc-gz "data/host_removed/${N}" > "data/bowtie/${N}_mapped_unmapped.sam"
done
