#!/bin/bash

# Usage: mapper.sh <read1.fastq.gz> <read2.fastq.gz> <output_prefix>
# This script aligns paired-end reads to the PhiX genome using Bowtie2. It is invoked separately for each pair of
# FASTQ files. The output is a SAM file containing the alignment results. This file needs to go in an output directory
# outside of the vendor outputs. In the next step, the SAM files will be merged, converted to BAM, sorted, and indexed.

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PHIX_INDEX_DIR="$SCRIPT_DIR/../../data"
PHIX_INDEX="$PHIX_INDEX_DIR/phix_index"

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <read1.fastq.gz> <read2.fastq.gz> <output_prefix>"
    exit 1
fi

READ1=$1
READ2=$2
OUTPUT_PREFIX=$3

if [ ! -f "${PHIX_INDEX}.1.bt2" ]; then
    echo "Error: PhiX index files not found in $PHIX_INDEX_DIR"
    exit 1
fi

echo "Aligning reads to PhiX genome..."
bowtie2 -x "$PHIX_INDEX" -1 "$READ1" -2 "$READ2" -S "${OUTPUT_PREFIX}.sam"

echo "Mapping complete for $OUTPUT_PREFIX"