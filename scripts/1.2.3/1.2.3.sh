# Metric: Assess sequencing quality using spiked-in PhiX control

# Align reads to PhiX genome
# Align PhiX reads:
# First, you need to identify and align the PhiX reads in your dataset.
# You can do this using a tool like Bowtie2 or BWA.
bowtie2 -x phix_index -U your_fastq_file.fq -S phix_aligned.sam

# Calculate alignment statistics:
# Use samtools to convert the SAM file to BAM, sort it, and generate alignment statistics.
samtools view -bS phix_aligned.sam | samtools sort - -o phix_aligned_sorted.bam
samtools flagstat phix_aligned_sorted.bam > phix_alignment_stats.txt

# Calculate error rate and other metrics:
# You can use a tool like Picard's CollectAlignmentSummaryMetrics to get detailed metrics:
java -jar picard.jar CollectAlignmentSummaryMetrics R=phix_genome.fasta I=phix_aligned_sorted.bam O=phix_metrics.txt

# Targeted analysis
# * a) PhiX Error Rate (2 points)
# * b) PhiX Alignment Rate (1 point)
# * c) PhiX Coverage Uniformity (1 point)
python assess_phix_quality.py phix_aligned_sorted.bam