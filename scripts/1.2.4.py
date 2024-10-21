import gzip
from collections import defaultdict
import sys


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

            # Print progress every million reads
            if total_reads % 1000000 == 0:
                print(f"Processed {total_reads} reads...", file=sys.stderr)

    # Calculate duplicate reads and rate
    duplicate_reads = sum(count - 1 for count in read_counts.values() if count > 1)
    duplicate_rate = duplicate_reads / total_reads if total_reads > 0 else 0

    print(f"Total reads: {total_reads}")
    print(f"Unique reads: {len(read_counts)}")
    print(f"Duplicate reads: {duplicate_reads}")
    print(f"Duplicate Read Rate: {duplicate_rate:.4f}")
    print(f"Percentage of duplicate reads: {duplicate_rate * 100:.2f}%")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.fastq[.gz]")
        sys.exit(1)

    calculate_duplicate_rate(sys.argv[1])