#!/usr/bin/env python

import pysam
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

def calculate_gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq)

def analyze_bam(bam_file, reference_file, window_size=100, step_size=50):
    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Load reference sequence
    reference = list(SeqIO.parse(reference_file, "fasta"))[0]

    # Initialize results list
    results = []

    # Iterate over the reference sequence in windows
    for i in range(0, len(reference.seq) - window_size + 1, step_size):
        window_seq = reference.seq[i:i+window_size]
        gc_content = calculate_gc_content(window_seq)

        # Calculate coverage for this window
        coverage = bam.count_coverage(reference.id, i, i+window_size)
        avg_coverage = np.mean([np.mean(cov) for cov in coverage])

        results.append((i, i+window_size, gc_content, avg_coverage))

    bam.close()

    return results

def main(bam_file, reference_file, output_file):
    results = analyze_bam(bam_file, reference_file)

    # Write results to CSV
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Window Start', 'Window End', 'GC Content', 'Average Coverage'])
        writer.writerows(results)

    # Calculate bias
    gc_rich_coverage = np.mean([r[3] for r in results if r[2] > 0.5])
    at_rich_coverage = np.mean([r[3] for r in results if r[2] <= 0.5])
    bias = gc_rich_coverage / at_rich_coverage if at_rich_coverage > 0 else float('inf')

    print(f"GC-rich regions average coverage: {gc_rich_coverage:.2f}")
    print(f"AT-rich regions average coverage: {at_rich_coverage:.2f}")
    print(f"Coverage bias (GC-rich / AT-rich): {bias:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <bam_file> <reference_fasta> <output_csv>")
        sys.exit(1)

    bam_file = sys.argv[1]
    reference_file = sys.argv[2]
    output_file = sys.argv[3]

    main(bam_file, reference_file, output_file)