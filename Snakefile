import pandas as pd
shell.executable("bash")

"""
Workflow is divided into three steps:
1. Trim adapter and low quality bases from short sequence reads (FASTQ format)
2. Align reads to a reference genome
3. Merge reads originating from the same biological sample
"""

# Config file
configfile: "config.yaml"

# Units file is a table that maps unit to sample information
# Here, a unit is any independent combination of biological sample, library prep, and sequencing run
# For variant calling from DNAseq data, it is often desirable to QC at the level of units
# After QC, units originating from the biological samples are merged

# units.tsv maps unit to sample info and specifies single- vs. paired-end reads
# reference for this approach: https://github.com/snakemake-workflows/rna-seq-star-deseq2
units = pd.read_csv(config["units"], index_col=["sample","unit"], dtype=str, sep = "\t")
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

# To merge sample info make a simple dictionary that maps sample to unit info, with the
# path specified. This isn't necessary and could easily be coded more elegantly.
# But it works.
samples = {}

# Here, merge after alignment. With more QC steps, the path at which you execute the merge
# will likely change
for i in units.index.levels[0]:
    samples[i] = ["data/aligned/{}-{}.bam".format(i,j) for j in units.loc[i].index]

# function determines which units are single-end reads
def is_single_end(sample,unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

# function returns reads as input for trim step with snakemake wildcards for sample and unit
def get_fastq(wildcards):
    return "data/reads/"+units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

# function returns reads as input for alignment step with snakemake wildcards for sample and unit
def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

# master rule
rule all:
    input:
        ["data/merged/{}.bam".format(x) for x in units.index.levels[0]],
        config["genome"]+".bwt"

# trim adapter and low quality bases from single end reads
rule cutadapt:
    input:
        get_fastq
    output:
        fastq="data/trimmed/{sample}-{unit}.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    threads: config["cutadapt"]["threads"]
    params:
        "a {} -q {}".format(config["cutadapt"]["adapter"], config["cutadapt"]["qual"])
    log:
        "log/trim/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/se"

# trim adapter and low quality bases from paired end reads
rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="data/trimmed/{sample}-{unit}-1.fastq.gz",
        fastq2="data/trimmed/{sample}-{unit}-2.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    threads:
        config["cutadapt"]["threads"]
    params:
        "-a {} -A {} -q {}".format(config["cutadapt"]["adapter"], config["cutadapt"]["adapter"], config["cutadapt"]["qual"])
    log:
        "log/trim/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"

# index reference genome
rule ref_index:
    input:
        config["genome"]
    output:
        config["genome"]+".bwt"
    log:
        "log/ref_index/index.log"
    shell:
        "bwa index {input}"

# align reads to reference genome
# here, alignments are piped directly to samtools sort
# 3:1 thread balancing might not be ideal. I've only halfheartedly tested this.
rule align:
    input:
        reads=get_trimmed,
        ref=config["genome"],
        index=config["genome"]+'.bwt'
    output:
        "data/aligned/{sample}-{unit}.bam"
    log:
        "log/align/{sample}-{unit}.log"
    threads:
        config["align"]["threads"]
    params:
        rg="'@RG\\t:ID:{unit}\\tSM:{sample}'",
        bwa_threads=3*config["align"]["threads"] // 4,
        sort_threads=config["align"]["threads"] // 4,
        sort_mem=config["align"]["sort_mem"]
    shell: """
        bwa mem -R {params.rg} -t {params.bwa_threads} {input.ref} {input.reads} 2> {log} | \
        samtools sort -@{params.sort_threads} -m {params.sort_mem} -o {output} -
    """

# merge alignments by biological sample
rule merge:
    input:
        lambda x: samples[x.sample]
    output:
        "data/merged/{sample}.bam"
    shell:
        "samtools merge {output} {input}"