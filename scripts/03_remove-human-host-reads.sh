#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N remove_host
#PBS -e remove_host.err
#PBS -o remove_host.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=48:00:00

# Load modules
module load tools
module load bowtie2/2.5.2

# Setup working directory and paths
BASE="path to working directory"
mkdir -p $BASE/results/{host_removed,bowtie}
samples=$(cat "Path to sample list file")

# Build Bowtie2 index for the human reference genome (GRCh38)
bowtie2-build $BASE/data/db/ncbi_dataset/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna ${path}/data/db/GRCh38_noalt_as

# Map trimmed reads to human reference (GRCh38) and extract unmapped (host-removed) reads
for sample in $samples; do
    f="$BASE/results/trimmomatic/${sample}_trimmed_1.fq.gz"; 
    r="$BASE/results/trimmomatic/${sample}_trimmed_2.fq.gz"
    echo 
    echo 'processing' $sample
    bowtie2 -p 40 \
            -x $BASE/data/db/GRCh38_noalt_as \
            -1 $f -2 $r \
            --un-conc-gz "$BASE/results/host_removed/${sample}_%.fq.gz" > "$BASE/results/bowtie/${sample}_mapped_unmapped.sam"
done
