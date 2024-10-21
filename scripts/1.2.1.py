#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO
from collections import defaultdict


def calculate_q_score_decay(seq, qual):
    homopolymer_runs = defaultdict(list)
    current_base = seq[0]
    run_start = 0

    for i, base in enumerate(seq):
        if base != current_base:
            if i - run_start >= 4:  # Only consider runs of 4 or more
                homopolymer_runs[current_base].append((run_start, i))
            current_base = base
            run_start = i

    # Check for a run that ends at the sequence end
    if len(seq) - run_start >= 4:
        homopolymer_runs[current_base].append((run_start, len(seq)))

    total_decay = 0
    total_comparisons = 0

    for runs in homopolymer_runs.values():
        for start, end in runs:
            q_scores = [ord(q) - 33 for q in qual[start:end]]  # Convert to numeric Q scores
            for i in range(1, len(q_scores)):
                decay = q_scores[i - 1] - q_scores[i]
                total_decay += decay
                total_comparisons += 1

    return total_decay / total_comparisons if total_comparisons > 0 else 0


def process_fastq(file_path):
    total_decay = 0
    total_reads = 0

    with gzip.open(file_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            decay = calculate_q_score_decay(str(record.seq), record.letter_annotations["phred_quality"])
            total_decay += decay
            total_reads += 1

            # Print progress every 10000 reads
            if total_reads % 10000 == 0:
                print(f"Processed {total_reads} reads...", file=sys.stderr)

    return total_decay / total_reads if total_reads > 0 else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.fastq")
        sys.exit(1)

    avg_decay = process_fastq(sys.argv[1])
    print(f"Average Q score decay per additional base in homopolymer runs: {avg_decay:.4f}")