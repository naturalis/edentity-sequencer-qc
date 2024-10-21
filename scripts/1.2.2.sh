#!/bin/bash

# Check if a file was provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input.fastq.gz>"
    exit 1
fi

input_file="$1"

# Check if the file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# Count the number of reads (divide by 4 as each FASTQ record has 4 lines)
nreads=$(zcat "$input_file" | awk 'END {print NR/4}')

# Get the base filename without path
filename=$(basename "$input_file")

# Determine the category based on the filename
if [[ "$filename" == *"Empty"* ]]; then
    category="Empty"
elif [[ "$filename" == *"Undetermined"* ]]; then
    category="Undetermined"
else
    category="Demultiplexed"
fi

# Output the result
echo "$category,$nreads"