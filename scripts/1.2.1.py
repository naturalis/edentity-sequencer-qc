#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO
from collections import defaultdict


def calculate_q_score_decay(seq, qual):
    """
    Calculate the average Q score decay per additional base in homopolymer runs.
    :param seq: a Seq record
    :param qual: a quality string
    :return:
    """
    homopolymer_runs = defaultdict(list)
    current_base = seq[0]
    run_start = 0

    # iterates over the bases in the sequence, keeping track of their position i
    for i, base in enumerate(seq):

        # there is a change relative to the focal base
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

    # iterate over all homopolymer runs
    for runs in homopolymer_runs.values():
        for start, end in runs:
            q_scores = qual[start:end]  # Quality scores are already numeric
            for i in range(1, len(q_scores)):
                decay = q_scores[i - 1] - q_scores[i]
                total_decay += decay
                total_comparisons += 1

    return total_decay / total_comparisons if total_comparisons > 0 else 0


def process_fastq(file_path):
    """
    Calculate the average Q score decay per additional base in homopolymer runs for a gzipped FASTQ file.
    :param file_path: a gzipped FASTQ file
    :return:
    """
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
        print("Usage: python script.py input.fastq.gz")
        sys.exit(1)

    avg_decay = process_fastq(sys.argv[1])
    print(f"Average Q score decay per additional base in homopolymer runs: {avg_decay:.4f}")