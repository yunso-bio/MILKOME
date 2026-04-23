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

# Setup working directory and paths
path="path to working directory"
mkdir -p $path/results/bins/avamb/{catalogue,bams,sorted,bins}
stool_samples=$(cat "path to a list of stool samples")

# Run Avamb for stool samples with vae-aaee model
for N in $stool_samples; do
    assembly="$path/result/megahit/${N}/final.contigs.fa"
    f="$path/results/host_removed/${N}_1.fq.gz"
    r="$path/results/host_removed/${N}_2.fq.gz"
    bam="$path/results/bins/avamb/bams/${N}.bam"
    sorted="$path/results/bins/avamb/sorted/${N}_sorted.bam"
    mmi="$path/results/bins/avamb/catalogue/${N}.mmi"

    # Index sample
    echo 'processing' $N 
    minimap2 -t 40 -d $mmi $assembly
    minimap2 -t 40 -N 5 -ax sr $mmi --split-prefix mmsplit $f $r | samtools view -F 3584 -b --threads 40 > $bam
    echo 'processing' $N
    samtools sort $bam -@ 40 -o $sorted 
    
    # Run Vamb
    vamb --outdir $path/result/bins/avamb/bins/${N} --model vae-aae --fasta $assembly --bamfiles $sorted --minfasta 200000
done

# Run Avamb for enriched samples with vae-aaee model
enriched_samples=$(cat "path to a list of stool samples")

for N in $enriched_samples; do
    assembly="$path/result/megahit/${N}/final.contigs.fa"
    f="$path/results/trimmomatic/${N}_trimmed_1.fq.gz" 
    r="$path/results/trimmomatic/${N}_trimmed_2.fq.gz"
    bam="$path/results/bins/avamb/bams/${N}.bam"
    sorted="$path/results/bins/avamb/sorted/${N}_sorted.bam"
    mmi="$path/results/bins/avamb/catalogue/${N}.mmi"

    # Index sample
    echo 'processing' $N 
    minimap2 -t 40 -d $mmi $assembly
    minimap2 -t 40 -N 5 -ax sr $mmi --split-prefix mmsplit $f $r | samtools view -F 3584 -b --threads 40 > $bam
    echo 'processing' $N
    samtools sort $bam -@ 40 -o $sorted 
    
    # Run Vamb
    vamb --outdir $path/results/bins/avamb/bins/${N} --model vae-aae --fasta $assembly --bamfiles $sorted --minfasta 200000
done
