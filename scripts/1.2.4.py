#!/usr/bin/env python
import gzip
from collections import defaultdict
import sys
import logging
import os

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', stream=sys.stderr)

def calculate_duplicate_rate(fastq_file):
    read_counts = defaultdict(int)
    total_reads = 0

    # Determine if the file is gzipped
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'

    with open_func(fastq_file, mode) as f:
        while True:
            # Read the four lines of each FASTQ record
            identifier = f.readline().strip()
            if not identifier:  # End of file
                break
            sequence = f.readline().strip()
            f.readline()  # Skip the '+' line
            f.readline()  # Skip the quality line

            # Count the occurrences of each unique sequence
            read_counts[sequence] += 1
            total_reads += 1

            # Log progress every million reads
            if total_reads % 1000000 == 0:
                logging.info(f"Processed {total_reads} reads...")

    # Calculate duplicate reads and rate
    duplicate_reads = sum(count - 1 for count in read_counts.values() if count > 1)
    duplicate_rate = duplicate_reads / total_reads if total_reads > 0 else 0

    # Prepare CSV output
    file_name = os.path.basename(fastq_file)
    unique_reads = len(read_counts)
    duplicate_percentage = duplicate_rate * 100

    # Print CSV output
    print(f"{file_name},{total_reads},{unique_reads},{duplicate_reads},{duplicate_rate:.4f},{duplicate_percentage:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.fastq[.gz]")
        sys.exit(1)

    calculate_duplicate_rate(sys.argv[1])