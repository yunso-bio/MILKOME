#!/bin/bash

# Check if the file is provided as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <vamb filename>"
    exit 1
fi

# Get the filename from the argument
filename=$1

# Check if the file exists
if [ ! -f "$filename" ]; then
    echo "File not found!"
    exit 1
fi

# Get the directory of the input file
directory=$(dirname "$filename")

# Generate a new file with swapped columns in the same directory
awk -F'\t' '{ temp = $1; $1 = $2; $2 = temp; print }' OFS='\t' "$filename" > "$directory/swapped_$(basename "$filename")"

echo "Swapped columns and saved to $directory/swapped_$(basename "$filename")"
