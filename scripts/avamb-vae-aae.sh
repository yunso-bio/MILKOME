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

mkdir -p bins/avamb/{catalogue,bams,sorted,bins}

samples="a list of sampleIDs"

# Run Avamb with vae-aaee model
for N in $samples; do
    assembly="data/megahit/${N}/final.contigs.fa"
    f="data/trimmomatic/${N}_trimmed_1.fq.gz" 
    r="data/trimmomatic/${N}_trimmed_2.fq.gz"
    bam="bins/avamb/bams/${N}.bam"
    sorted="bins/avamb/sorted/${N}_sorted.bam"
    mmi="bins/avamb/catalogue/${N}.mmi"

    # Index sample
    echo 'processing' $N 
    minimap2 -t 40 -d $mmi $assembly
    minimap2 -t 40 -N 5 -ax sr $mmi --split-prefix mmsplit $f $r | samtools view -F 3584 -b --threads 40 > $bam
    echo 'processing' $N
    samtools sort $bam -@ 40 -o $sorted 
    
    # Run Vamb
    vamb --outdir bins/avamb/bins/${N} --model vae-aae --fasta $assembly --bamfiles $sorted --minfasta 200000
done


mail -s 'AVAMB - binning is done' yunso@dtu.dk
