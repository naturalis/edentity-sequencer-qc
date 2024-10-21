#!/bin/bash

# Metric: Assess sequencing quality using spiked-in PhiX control

# Get the directory of the script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Set the path to the PhiX index relative to the script location
PHIX_INDEX_DIR="$SCRIPT_DIR/../../data"
PHIX_INDEX="$PHIX_INDEX_DIR/phix_index"

# Check if input files are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <read1.fastq.gz> <read2.fastq.gz>"
    exit 1
fi

READ1=$1
READ2=$2

# Check if PhiX index files exist
if [ ! -f "${PHIX_INDEX}.1.bt2" ]; then
    echo "Error: PhiX index files not found in $PHIX_INDEX_DIR"
    echo "Please ensure the PhiX index files (phix_index.*) are present in the data directory."
    exit 1
fi

# Align paired-end reads to PhiX genome
echo "Aligning reads to PhiX genome..."
bowtie2 -x "$PHIX_INDEX" -1 "$READ1" -2 "$READ2" -S phix_aligned.sam

# Calculate alignment statistics
echo "Calculating alignment statistics..."
samtools view -bS phix_aligned.sam | samtools sort - -o phix_aligned_sorted.bam
samtools index phix_aligned_sorted.bam
samtools flagstat phix_aligned_sorted.bam > phix_alignment_stats.txt

# Generate coverage statistics
echo "Generating coverage statistics..."
samtools coverage phix_aligned_sorted.bam > phix_coverage.txt

# Targeted analysis
echo "Performing targeted analysis..."
# * a) PhiX Error Rate (2 points)
# * b) PhiX Alignment Rate (1 point)
# * c) PhiX Coverage Uniformity (1 point)
python "$SCRIPT_DIR/assess_phix_quality.py" phix_aligned_sorted.bam phix_coverage.txt

# Clean up intermediate files
echo "Cleaning up..."
rm phix_aligned.sam