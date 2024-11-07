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

# Get the directory of the current script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Path to the mapper.sh script
MAPPER_SCRIPT="$SCRIPT_DIR/mapper.sh"

# Check if mapper.sh exists and is executable
if [ ! -x "$MAPPER_SCRIPT" ]; then
    echo "Error: mapper.sh not found or not executable in $SCRIPT_DIR"
    exit 1
fi

# Find and process all R1 files, then find their R2 counterparts
find "$INPUT_DIR" -name "*_R1.fastq.gz" | sort | while read r1_file; do
    r2_file="${r1_file/_R1/_R2}"

    # Check if R2 file exists
    if [ ! -f "$r2_file" ]; then
        echo "Warning: No matching R2 file found for $r1_file"
        continue
    fi

    # Extract the prefix for the output SAM file
    filename=$(basename "$r1_file")
    prefix="${filename/_R1.fastq.gz/}"
    output_prefix="$OUTPUT_DIR/$prefix"

    echo "Processing: $r1_file and $r2_file"
    echo "Output prefix: $output_prefix"

    # Call the external mapper.sh script
    "$MAPPER_SCRIPT" "$r1_file" "$r2_file" "$output_prefix"

    echo "Alignment complete for $prefix"
    echo "-----------------------------"
done

echo "All files processed."