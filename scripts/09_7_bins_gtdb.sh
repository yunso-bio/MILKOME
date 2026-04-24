#!/bin/bash

#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N gtdb
#PBS -e gtdb.err
#PBS -o gtdb.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=24:00:00

# Load modules
module load tools
module load gcc/12.2.0
module load fastani/1.34
module load dashing/0.4.0
module load gtdbtk/2.3.2

base_dir="path to a working directory"

# Run gtdbtk
gtdbtk classify_wf --cpus 40 \
				   --genome_dir $base_dir/bins/ani98 \
				   -x "fna" \
				   --out_dir  $base_dir/bins/gtdb \
				   --skip_ani_screen;
