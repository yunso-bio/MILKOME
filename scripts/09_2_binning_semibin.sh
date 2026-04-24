#!/bin/bash
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N easy-bin
#PBS -e easy-bin.err
#PBS -o easy-bin.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=48:00:00

set -euo pipefail

# Load modules
module load tools
module load semibin/2.1.0
module load samtools/1.18
module load minimap2/2.24r1122
module load bowtie2/2.5.2

# setup working direcoty and paths
BASE="path to working directory"
mkdir -p "$BASE/results/bins/semi_single_bin"/{maps,sams,bams,sorted,bins}

### Mapping ###
mapfile -t samples < "path to a list of samples"

for N in "${samples[@]}"; do
    echo "Processing $N"

    contig="$BASE/results/bins/semi_single_bin/maps/${N}_contigs"
    assembly="path/results/megahit/final.contigs.fa"
    f="$BASE/results/trimmomatic/${N}_trimmed_1.fq.gz"
    r="$BASE/results/trimmomatic/${N}_trimmed_2.fq.gz"

    sam="$BASE/results/bins/semi_single_bin/sams/${N}.sam"
    bam="$BASE/results/bins/semi_single_bin/bams/${N}.bam"
    mbam="$BASE/results/bins/semi_single_bin/bams/${N}.mapped.bam"
    sorted="$BASE/results/bins/semi_single_bin/sorted/${N}.mapped.sorted.bam"

    # Input checks
    if [[ ! -f "$assembly" || ! -f "$f" || ! -f "$r" ]]; then
        echo "Missing input for $N, skipping..."
        continue
    fi

    ### Index + mapping ###
    bowtie2-build -f "$assembly" "$contig"

    bowtie2 -q --fr \
        -x "$contig" \
        -1 "$f" \
        -2 "$r" \
        -S "$sam" \
        -p 40

    samtools view -h -b -S "$sam" -o "$bam" -@ 40
    samtools view -b -F 4 "$bam" -o "$mbam" -@ 40
    samtools sort "$mbam" -o "$sorted" -@ 40

    ### Binning ###
    outdir="$BASE/results/bins/semi_single_bin/bins/${N}"
    mkdir -p "$outdir"

    SemiBin2 single_easy_bin \
        -t 40 \
        --environment human_gut \
        -i "$assembly" \
        -b "$sorted" \
        -o "$outdir"
done
