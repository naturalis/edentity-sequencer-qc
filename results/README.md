# Results

This folder contains the intermediate results files. They are organized as follows:

## 1.1.1 Fraction of bases ≥ Q30 

This was computed using [1.1.1.sh](../scripts/1.1.1.sh), which produces CSV output with the following columns:

1. file name
2. total reads
3. average percentage of bases ≥ Q30 

* [Result for elements](elements-1.1.1-qc.csv)
* [Result for illumina](illumina-1.1.1-qc.csv)

## 1.1.2 Read Length Retention after Q30 Trimming

This was computed with [cutadapt.smk](../scripts/1.1.2/cutadapt.smk), which trims the reads using
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) at Q30. We assess the effect of this by comparing
FastQC reports before and after trimming.

## 1.2.1 Homopolymer Accuracy 

This was computed using [1.2.1.py](../scripts/1.2.1.py). This script assumes that a homopolymer
is a stretch of ≥4 identical bases, and produces CSV output with the following columns:

1. file name
2. total reads
3. average decay

## 1.2.2 Demultiplexing Efficiency 

This was computed using [1.2.2.sh](../scripts/1.2.2.sh), which produces CSV output with the following columns:

1. Category, i.e. Empty, Undetermined, or Demultiplexed
2. Number of reads

(Every line in the CSV is an input file.)

* [Result for elements](elements-1.2.2-demux.csv)
* [Result for illumina](illumina-1.2.2-demux.csv)

## 1.2.3 PhiX Control Performance 

This was performed by mapping all paired file sets to the PhiX reference genome in the [data](../data) folder
using bowtie2. This yielded a SAM file for each file pair. All SAM files are then merged, sorted and indexed
using samtools. Subsequently, PhiX Error Rate (2 points), PhiX Alignment Rate (1 point) and PhiX Coverage Uniformity 
(1 point) were computed with [assess_phix_quality.py](../scripts/1.2.3/assess_phix_quality.py).

## 1.2.4 Duplicate Read Rate 

This was computed using [1.2.4.py](../scripts/1.2.4.py), which produces CSV output with the following columns:

1. file name
2. total reads
3. unique reads
4. duplicate reads
5. duplicate rate
6. duplicate percentage

* [Result for elements](elements-1.2.4-dedup.csv)