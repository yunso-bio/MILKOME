#!/bin/bash
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N avamb
#PBS -e avamb.err
#PBS -o avamb.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=72:00:00

#Load modules 
module load tools
module load minimap2/2.24r1122
module load samtools/1.18
module load vamb/4.1.3 
module load checkm2/1.0.2

# Setup working directory and paths
base_dir="path to working directory"
mkdir -p $base_dir/results/bins/avamb/{catalogue,bams,sorted,bins,checkm2}

stool_samples=$(cat "path to a list of stool samples")

########################################
# STOOL SAMPLES
########################################
for N in $stool_samples; do
    assembly="$base_dir/result/megahit/${N}/final.contigs.fa"
    f="$base_dir/results/host_removed/${N}_1.fq.gz"
    r="$base_dir/results/host_removed/${N}_2.fq.gz"
    bam="$base_dir/results/bins/avamb/bams/${N}.bam"
    sorted="$base_dir/results/bins/avamb/sorted/${N}_sorted.bam"
    mmi="$base_dir/results/bins/avamb/catalogue/${N}.mmi"
    outdir="$base_dir/results/bins/avamb/bins/${N}"
    checkm_out="$base_dir/results/bins/avamb/checkm2/${N}"

    echo "Processing $N"

    # Index + mapping
    minimap2 -t 40 -d $mmi $assembly
    minimap2 -t 40 -N 5 -ax sr $mmi --split-prefix mmsplit $f $r | \
        samtools view -F 3584 -b --threads 40 > $bam

    samtools sort $bam -@ 40 -o $sorted 

    # Run VAMB
    vamb --outdir $outdir \
         --model vae-aae \
         --fasta $assembly \
         --bamfiles $sorted \
         --minfasta 200000

    # Run CheckM2 on bins
    echo "Running CheckM2 for $N"
    checkm2 predict \
        --input $outdir/bins \
        --output-directory $checkm_out \
        --threads 40 \
        --database_dir_path "$base_dir_dir/db/CheckM2_db/uniref100.KO.1.dmnd"
done

########################################
# ENRICHED SAMPLES
########################################
enriched_samples=$(cat "path to a list of enriched samples")

for N in $enriched_samples; do
    assembly="$base_dir/result/megahit/${N}/final.contigs.fa"
    f="$base_dir/results/trimmomatic/${N}_trimmed_1.fq.gz" 
    r="$base_dir/results/trimmomatic/${N}_trimmed_2.fq.gz"
    bam="$base_dir/results/bins/avamb/bams/${N}.bam"
    sorted="$base_dir/results/bins/avamb/sorted/${N}_sorted.bam"
    mmi="$base_dir/results/bins/avamb/catalogue/${N}.mmi"
    outdir="$base_dir/results/bins/avamb/bins/${N}"
    checkm_out="$base_dir/results/bins/avamb/checkm2/${N}"

    echo "Processing $N"

    # Index + mapping
    minimap2 -t 40 -d $mmi $assembly
    minimap2 -t 40 -N 5 -ax sr $mmi --split-prefix mmsplit $f $r | \
        samtools view -F 3584 -b --threads 40 > $bam

    samtools sort $bam -@ 40 -o $sorted 

    # Run VAMB
    vamb --outdir $outdir \
         --model vae-aae \
         --fasta $assembly \
         --bamfiles $sorted \
         --minfasta 200000

    # Run CheckM2
    echo "Running CheckM2 for $N"
    checkm2 predict \
        --input $outdir/bins \
        --output-directory $checkm_out \
        --threads 40 \
        --database_dir_path "$base_dir_dir/db/CheckM2_db/uniref100.KO.1.dmnd"
done
