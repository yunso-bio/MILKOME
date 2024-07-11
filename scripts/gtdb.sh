#!/bin/bash

#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N gtdb
#PBS -e gtdb.err
#PBS -o gtdb.log
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=24:00:00

module load tools
module load gcc/12.2.0
module load fastani/1.34
module load dashing/0.4.0
module load gtdbtk/2.3.2

cd /home/projects/dtu_00032/analysis/milk-cohort/bins/new-ani98

 
gtdbtk classify_wf --cpus 40 \
				   --genome_dir bins \
				   -x "fna" \
				   --out_dir gtdb \
				   --skip_ani_screen;

mail -s 'GTDB is done' yunso@dtu.dk
