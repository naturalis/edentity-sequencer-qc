
import os
import re
import glob

configfile: "/data/luka.lenaroto/edentity-sequencer-qc/configs/fastqc.yaml"
# constract the base names for the wildcards

SAMPLES,  EXTENSION = glob_wildcards(os.path.join(config["raw_data"], '{sample}_{extension}'))
# take unique parts of the elements
SAMPLES = list(set(SAMPLES))
#print(SAMPLES)

#remove the extension
SAMPLE_NAMES = []
skip = ["mappings", "fastq-before-trimming"]

for samp in SAMPLES:
    if samp.endswith("_R1") or samp.endswith("_R2"):
    
        SAMPLE_NAMES.append("_".join(samp.split("_")[:-1]))


print(SAMPLE_NAMES)
#forward_extension, reverse_extension = sorted(list(set(EXTENSION)), reverse=False)
rule all:
    input:
        expand(os.path.join(config["output-dir"], f"{{sample}}_R1_001_fastqc.html"), sample=SAMPLE_NAMES),
        expand(os.path.join(config["output-dir"], f"{{sample}}_R2_001_fastqc.zip"), sample=SAMPLE_NAMES)


rule fastqc:
    input:
        r1 = os.path.join(config['raw_data'], f"{{sample}}_R1_001.fastq.gz"),
        r2 = os.path.join(config['raw_data'],f"{{sample}}_R2_001.fastq.gz")
        
    output:
        report = os.path.join(config["output-dir"], f"{{sample}}_R1_001_fastqc.html"),
        zip_file = os.path.join(config["output-dir"], f"{{sample}}_R2_001_fastqc.zip")

    
    conda: "/data/luka.lenaroto/edentity-sequencer-qc/environment.yml"
    #benchmark: os.path.join(config["output-dir"], "cutadapt.log")

    params: config["output-dir"]


    shell:
        """
        fastqc {input.r1} {input.r2} -o {params}
        """
        
#