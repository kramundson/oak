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
	["data/depths/{}_Q20_depth.bed".format(x) for x in units.index.levels[0]],
        ["data/summarized_depth/{}_depth_summary.tsv".format(x) for x in units.index.levels[0]],
        config["genome"]+".bwt",
        config["genome"].split(".fa")[0] + "_GCN.tsv"

# create temporary folder for parallel-fastq-dump to use instead of /tmp/
rule tmp_folder:
    output:
        temp(touch("data/reads/tmp/flag"))

# run parallel fastq dump if sra files present but not reads
rule fastq_dump:
    input:
        "{folder}/{{id}}.sra".format(folder=config["fastq_dump"]["sra_location"].rstrip("/")),
        "data/reads/tmp/flag"
    output:
        "data/reads/{id}_pass_1.fastq.gz",
        "data/reads/{id}_pass_2.fastq.gz"
    params:
        # alternative temporary directory
        tmpdir="data/reads/tmp/",
        # native fastq-dump options
        options="--outdir data/reads/ --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip"
    threads:
        config["fastq_dump"]["threads"]
    log:
        "log/fastq_dump/{id}.log"
    shell:
        # --readids option will break bwa mem
        "parallel-fastq-dump --sra-id {input[0]} -t {threads} --tmpdir {params.tmpdir} {params.options} > {log} 2>&1"

# trim adapter and low quality bases from single end reads
rule cutadapt:
    input:
        get_fastq
    output:
        fastq="data/trimed/{sample}-{unit}.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    threads: config["cutadapt"]["threads"]
    params:
        "-a {} -q {}".format(config["cutadapt"]["adapter"], config["cutadapt"]["qual"])
    log:
        "log/trimmed/{sample}-{unit}.log"
    shell: """
        cutadapt {params} -j {threads} -o {output.fastq} {input} > {output.qc} 2> {log}
    """

# trim adapter and low quality bases from paired end reads
rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="data/trimmed/{sample}-{unit}-1.fastq.gz",
        fastq2="data/trimmed/{sample}-{unit}-2.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} -A {} -q {}".format(config["cutadapt"]["adapter"], config["cutadapt"]["adapter"], config["cutadapt"]["qual"])
    threads: config["cutadapt"]["threads"]
    log:
        "log/trimmed/{sample}-{unit}.log"
    shell: """
        cutadapt {params} -o {output.fastq1} -p {output.fastq2} -j {threads} {input} > {output.qc} 2> {log}
    """

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
        rg="'@RG\\tID:{unit}\\tSM:{sample}'",
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

rule depth:
    input:
        "data/merged/{sample}.bam"
    output:
        "data/depths/{sample}_Q20_depth.bed"
    log: "log/depth/{sample}.log"
    shell: """
        samtools depth -a -Q 20 {input} > {output} 2> {log}
    """

# Handling genome in top-level or in subfolders:

# Prepend "./" to ensure there will always be a slash, then rsplit on slash.
# Unpack tuple to place folder in first blank and the rest in last blank.
# If genome is in top-level, the folder will just be a dot.

# However relative file paths make snakemake mad
# so after populating the string, use split to strip the leading "./"

rule window:
    input:
        config["genome"]
    output:
        fai="{}.fai".format(config["genome"]),
        gen="{}.contigs".format(config["genome"]),
        win="{}/10k_{}.bed".format(*("./" + config["genome"]).rsplit("/", 1)).split("./", 1)[-1]
    shell: """
        samtools faidx {input}
        awk -v OFS='\t' '{{print $1,$2}}' {output.fai} > {output.gen}
        bedtools makewindows -w 10000 -g {output.gen} > {output.win}
    """

rule window_depth:
    input:
        win="{}/10k_{}.bed".format(*("./" + config["genome"]).rsplit("/", 1)).split("./", 1)[-1],
        dep="data/depths/{sample}_Q20_depth.bed"
    output:
        "data/summarized_depth/{sample}_depth_summary.tsv"
    log:
        "log/summarized_depth/{sample}.log"
    shell: """
        python scripts/window_depth_summarizer.py -w {input.win} -b {input.dep} -o {output} > {log} 2>&1
    """

rule GCN_content:
    input:
        win="{}/10k_{}.bed".format(*("./" + config["genome"]).rsplit("/", 1)).split("./", 1)[-1],
        genome=config["genome"]
    output:
        config["genome"].split(".fa")[0] + "_GCN.tsv"
    shell: """
        python scripts/gcn_content_summarizer.py {input.win} {input.genome} > {output}
    """
