#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N stool-mp
#PBS -e stool-mp.err
#PBS -o stool-mp.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=48:00:00


module load tools
module load bowtie2/2.5.2

cd /home/projects/dtu_00032/analysis/milk-cohort
#mkdir data/host_removed
#mkdir data/bowtie

file=$(cat logs/stools.txt)
for N in $file; do
    f="data/trimmomatic/${N}_trimmed_1.fq.gz"; 
    r="data/trimmomatic/${N}_trimmed_2.fq.gz"
    echo 
    echo 'processing' $N
    bowtie2 -p 40 \
            -x /home/projects/dtu_00032/db/GRCh38_noalt_as/GRCh38_noalt_as \
            -1 $f -2 $r \
            --un-conc-gz "data/host_removed/${N}" > "data/bowtie/${N}_mapped_unmapped.sam"
done


#bowtie2-build ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38
