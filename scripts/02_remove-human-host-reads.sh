#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N stool-mp
#PBS -e stool-mp.err
#PBS -o stool-mp.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=48:00:00

# Load modules
module load tools
module load bowtie2/2.5.2

# Setup working directory and paths
path="path to working directory"
cd "$path" || exit
mkdir -p data/{host_removed,bowtie}
samples=$(cat "Path to sample list file")

# Build Bowtie2 index for the human reference genome (GRCh38)
bowtie2-build ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38

# Map trimmed reads to human reference (GRCh38) and extract unmapped (host-removed) reads
for sample in $samples; do
    f="data/trimmomatic/${N}_trimmed_1.fq.gz"; 
    r="data/trimmomatic/${N}_trimmed_2.fq.gz"
    echo 
    echo 'processing' $sample
    bowtie2 -p 40 \
            -x /home/projects/dtu_00032/db/GRCh38_noalt_as/GRCh38_noalt_as \
            -1 $f -2 $r \
            --un-conc-gz "data/host_removed/${sample}" > "data/bowtie/${sample}_mapped_unmapped.sam"
done
