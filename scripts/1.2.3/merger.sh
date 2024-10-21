#!/bin/bash

# Usage: merger.sh <sam_file1> [<sam_file2> ...]
# This script merges multiple SAM files, converts the output to BAM, sorts it, and indexes it. It also calculates
# alignment statistics and coverage statistics for the merged BAM file.

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <sam_file1> [<sam_file2> ...]"
    exit 1
fi

echo "Merging SAM files..."
samtools merge -O SAM merged.sam "$@"

echo "Converting to BAM, sorting, and indexing..."
samtools view -bS merged.sam | samtools sort - -o phix_aligned_sorted.bam
samtools index phix_aligned_sorted.bam

echo "Calculating alignment statistics..."
samtools flagstat phix_aligned_sorted.bam > phix_alignment_stats.txt

echo "Generating coverage statistics..."
samtools coverage phix_aligned_sorted.bam > phix_coverage.txt

echo "Cleaning up..."
rm merged.sam