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
module load samtools/1.18
module load minimap2/2.24r1122
module load bowtie2/2.5.2
module load checkm2/1.0.2

# setup working direcoty and paths
base_dir="path to working directory"
mkdir -p "$base_dir/results/bins/semi_single_bin"/{maps,sams,bams,sorted,bins,checkm2}

mapfile -t samples < "/path/to/sample_list.txt"

### Mapping + binning ###
for N in "${samples[@]}"; do
    echo "Processing $N"

    assembly="$base_dir/results/megahit/${N}/final.contigs.fa"
    contig="$base_dir/results/bins/semi_single_bin/maps/${N}_contigs"

    f="$base_dir/results/trimmomatic/${N}_trimmed_1.fq.gz"
    r="$base_dir/results/trimmomatic/${N}_trimmed_2.fq.gz"

    sam="$base_dir/results/bins/semi_single_bin/sams/${N}.sam"
    bam="$base_dir/results/bins/semi_single_bin/bams/${N}.bam"
    mbam="$base_dir/results/bins/semi_single_bin/bams/${N}.mapped.bam"
    sorted="$base_dir/results/bins/semi_single_bin/sorted/${N}.mapped.sorted.bam"

    if [[ ! -f "$assembly" || ! -f "$f" || ! -f "$r" ]]; then
        echo "Missing input for $N, skipping..."
        continue
    fi

    ### Index + mapping ###
    bowtie2-build "$assembly" "$contig"

    bowtie2 -q --fr \
        -x "$contig" \
        -1 "$f" -2 "$r" \
        -S "$sam" \
        -p 40

    samtools view -h -b "$sam" -o "$bam" -@ 40
    samtools view -b -F 4 "$bam" -o "$mbam" -@ 40
    samtools sort "$mbam" -o "$sorted" -@ 40

    rm -f "$sam" "$bam"

    ### Binning ###
    outdir="$base_dir/results/bins/semi_single_bin/bins/${N}"
    mkdir -p "$outdir"

    SemiBin2 single_easy_bin \
        -t 40 \
        --environment human_gut \
        -i "$assembly" \
        -b "$sorted" \
        -o "$outdir"
done

### CheckM2 ###
for N in "${samples[@]}"; do
    bins="$base_dir/results/bins/semi_single_bin/bins/${N}/output_bins"

    for f in "$bins"/*.fa.gz; do
        i=$(basename "$f" .fa.gz)

        echo "Running CheckM2 for $i"

        checkm2 predict \
            --threads 40 \
            --input "$f" \
            -x fa.gz \
            --output-directory "$base_dir/results/bins/semi_single_bin/checkm2/${N}_${i}" \
            --database_path "$base_dir/db/CheckM2_database/uniref100.KO.1.dmnd"
    done
done
