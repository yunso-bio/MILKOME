#!/bin/bash
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N bakta
#PBS -e bakta.err
#PBS -o bakta.log
#PBS -l nodes=1:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=72:00:00

# Load modules
module load tools
module load bakta/1.9.3

# Paths
base_dir="path to a working directory"
input_dir="$base_dir/bins/das_ani98"
BTdb="$base_dir/db/bakta_db"
BTpath="$base_dir/bins/das_ani98/bakta"

mkdir -p "$BTpath"

# Read samples safely
mapfile -t samples_array < <(ls "$input_dir")

declare -A sample_map

# Assign IDs
for i in "${!samples_array[@]}"; do
    old_name="${samples_array[$i]}"
    new_name="A$(printf "%04d" $((i+1)))"
    sample_map["$old_name"]="$new_name"
    echo "Old name: $old_name -> New name: $new_name"
done

# Mapping file
mapping_file="$base_dir/bins/das-bins/sample_mapping.txt"
mkdir -p "$(dirname "$mapping_file")"
: > "$mapping_file"

for old_name in "${samples_array[@]}"; do
    echo "$old_name -> ${sample_map[$old_name]}" >> "$mapping_file"
done

# Run Bakta
for old_name in "${samples_array[@]}"; do
    N="${sample_map[$old_name]}"
    echo
    echo "Running bakta for $old_name"
    echo

    bakta --threads 40 \
          --db "$BTdb" \
          --compliant \
          --verbose \
          --output "$BTpath/$N" \
          --locus-tag "$N" \
          --prefix "$N" \
          "$input_dir/${old_name}"
done
