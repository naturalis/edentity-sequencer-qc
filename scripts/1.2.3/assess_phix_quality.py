#!/usr/bin/env python
import pysam
import sys
import numpy as np
from collections import defaultdict


def assess_phix_quality(bam_file):
    total_reads = 0
    mapped_reads = 0
    total_bases = 0
    mismatched_bases = 0
    coverage = defaultdict(int)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        phix_length = sum(bam.lengths)  # Assuming PhiX is the only reference

        for read in bam.fetch():
            total_reads += 1
            if read.is_unmapped:
                continue

            mapped_reads += 1
            total_bases += read.query_length

            for qpos, refpos in read.get_aligned_pairs(with_seq=True):
                if qpos is None or refpos is None:
                    continue  # Skip insertions/deletions

                query_base = read.query_sequence[qpos].upper()
                ref_base = refpos[2].upper()

                coverage[refpos] += 1

                if query_base != ref_base:
                    mismatched_bases += 1

    # a) PhiX Error Rate (2 points)
    error_rate = mismatched_bases / total_bases if total_bases > 0 else 0

    # b) PhiX Alignment Rate (1 point)
    alignment_rate = mapped_reads / total_reads if total_reads > 0 else 0

    # c) PhiX Coverage Uniformity (1 point)
    coverage_values = list(coverage.values())
    mean_coverage = np.mean(coverage_values)
    std_coverage = np.std(coverage_values)
    cv_coverage = std_coverage / mean_coverage if mean_coverage > 0 else 0

    # Calculate uniformity as the percentage of bases within ±20% of the mean coverage
    within_range = sum(0.8 * mean_coverage <= cov <= 1.2 * mean_coverage for cov in coverage_values)
    uniformity = within_range / len(coverage_values) if coverage_values else 0

    print(f"PhiX Error Rate: {error_rate:.6f}")
    print(f"PhiX Alignment Rate: {alignment_rate:.4f}")
    print(f"PhiX Coverage Uniformity:")
    print(f"  - Coefficient of Variation: {cv_coverage:.4f}")
    print(f"  - Percentage within ±20% of mean coverage: {uniformity:.4f}")
    print(f"\nAdditional Information:")
    print(f"  - Total Reads: {total_reads}")
    print(f"  - Mapped Reads: {mapped_reads}")
    print(f"  - Mean Coverage: {mean_coverage:.2f}")
    print(f"  - Coverage Range: {min(coverage_values)} - {max(coverage_values)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py phix_aligned_sorted.bam")
        sys.exit(1)

    assess_phix_quality(sys.argv[1])