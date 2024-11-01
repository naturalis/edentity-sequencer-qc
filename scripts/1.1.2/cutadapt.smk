
import os
import re
import glob

configfile: "../config/cutadapt.yaml"
# constract the base names for the wildcards

SAMPLES,  EXTENSION = glob_wildcards(os.path.join(config["raw_data"], '{sample}_{extension}'))
# take unique parts of the elements
SAMPLES = list(set(SAMPLES))

#remove the extension
SAMPLE_NAMES = []
skip = ["mappings", "fastq-before-trimming", "data"]

for samp in SAMPLES:
    if samp.endswith("_R1") or samp.endswith("_R2"):
    
        SAMPLE_NAMES.append("_".join(samp.split("_")[:-1]))


rule all:
    input:
        expand(os.path.join(config["cutadapt-output"], f'{{sample}}_trimmed_R1_001.fastq.gz'), sample=SAMPLE_NAMES),
        expand(os.path.join(config["cutadapt-output"], f"{{sample}}_trimmed_R2_001.fastq.gz"), sample=SAMPLE_NAMES)


rule cutadapt:
    input:
        r1 = os.path.join(config['raw_data'], f"{{sample}}_R1_001.fastq.gz"),
        r2 = os.path.join(config['raw_data'],f"{{sample}}_R2_001.fastq.gz")
        
    output:
        trimmed_1 = os.path.join(config["cutadapt-output"], f"{{sample}}_trimmed_R1_001.fastq.gz"),
        trimmed_2 = os.path.join(config["cutadapt-output"], f"{{sample}}_trimmed_R2_001.fastq.gz")

    
    conda: "/data/luka.lenaroto/edentity-sequencer-qc/environment.yml"
    #benchmark: os.path.join(config["cutadapt-output"], "cutadapt.log")

    shell:
        """
        cutadapt -q 30 -o {output.trimmed_1} {input.r1} && cutadapt -q 30 -o {output.trimmed_2} {input.r2} 
        """
        
#