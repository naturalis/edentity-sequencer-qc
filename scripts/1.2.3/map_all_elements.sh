#!/bin/bash

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Path to the PhiX index
PHIX_INDEX="../../data/phix_index"

# Function to process a pair of FASTQ files
process_pair() {
    local r1="$1"
    local r2="$2"
    local prefix="$3"

    echo "Processing: $r1 and $r2"
    echo "Output prefix: $prefix"

    # Run bowtie2 alignment
    bowtie2 -x "$PHIX_INDEX" -1 "$r1" -2 "$r2" -S "${prefix}.sam"

    echo "Alignment complete for $prefix"
    echo "-----------------------------"
}

# Find and process all R1 files, then find their R2 counterparts
find "$INPUT_DIR" -name "*_R1_001.fastq.gz" | sort | while read r1_file; do
    r2_file="${r1_file/_R1_/_R2_}"

    # Check if R2 file exists
    if [ ! -f "$r2_file" ]; then
        echo "Warning: No matching R2 file found for $r1_file"
        continue
    fi

    # Extract the prefix for the output SAM file
    filename=$(basename "$r1_file")
    prefix="${filename/_R1_001.fastq.gz/}"
    output_prefix="$OUTPUT_DIR/$prefix"

    # Process the pair
    mapper.sh "$r1_file" "$r2_file" "$output_prefix"
done

echo "All files processed."