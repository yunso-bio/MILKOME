#!/bin/bash

# Load modules 
module load tools
module load bakta/1.9.3

# path
path="path to working directory"
BTdb="$path/db/bakta_db"
BTpath="$path/results/bakta/"
mkdir -p $BTpath

# Map ids
samples=$(cat "path to a list of stool samples")
IFS=' ' read -r -a samples_array <<< "$samples"

declare -A sample_map

for i in "${!samples_array[@]}"; do
    old_name="${samples_array[$i]}"
    new_name="META$(printf "%04d" $((i+1)))"
    sample_map["$old_name"]="$new_name"
    echo "Old name: $old_name -> New name: $new_name"
done

# Create a mapping file
mapping_file="path_to/sample_mapping.txt"
: > "$mapping_file"
for old_name in "${!sample_map[@]}"; do
    echo "$old_name -> ${sample_map[$old_name]}" >> "$mapping_file"
done

# Run Bakta
for old_name in "${!sample_map[@]}"; do
    N="${sample_map[$old_name]}"
    echo 
    echo "Running bakta for $old_name"
    echo
    bakta --threads 40 \
          --db $BTdb \
          --compliant \
          --meta \
          --verbose \
          --output $BTpath/$N \
          --locus-tag $N \
          --prefix $N \
          "$path/results/megahit/${old_name}/final.contigs.fa"
done
