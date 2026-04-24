#!/bin/bash
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N coverm
#PBS -e coverm.err
#PBS -o coverm.log
#PBS -M yunso@dtu.dk
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=72:00:00

base_dir="path to bins"

### Coverm Clustering Dereplicate bins

# Load modules
module load tools
module load gcc/12.2.0
module load fastani/1.34
module load dashing/0.4.0
module load coverm/0.7.0

coverm cluster -t 40  \
				--genome-fasta-directory $base_dir/bins/filtered_bins \
				--output-cluster-definition ani98-cluster-definition \
				--checkm2-quality-report "$base_dir/bins/checkm2/concat_checkm_quality.tsv" \
				--min-completeness 50 \
				--max-contamination 5 \
				--ani 98 \
				--output-representative-fasta-directory $base_dir/bins/das_ani98;
