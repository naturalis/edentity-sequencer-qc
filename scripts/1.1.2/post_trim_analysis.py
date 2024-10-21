import sys
from Bio import SeqIO


def calculate_average_length(fastq_file):
    total_length = 0
    read_count = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        total_length += len(record.seq)
        read_count += 1
    return total_length / read_count if read_count > 0 else 0


def main(input_file, trimmed_file):
    original_avg_length = calculate_average_length(input_file)
    trimmed_avg_length = calculate_average_length(trimmed_file)

    retention = (trimmed_avg_length / original_avg_length) * 100 if original_avg_length > 0 else 0

    print(f"Original average read length: {original_avg_length:.2f}")
    print(f"Trimmed average read length: {trimmed_avg_length:.2f}")
    print(f"Read Length Retention after Q30 Trimming: {retention:.2f}%")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.fastq trimmed_output.fastq")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])