import pandas as pd
import re
from Bio import SeqIO
shell.executable("bash")

"""
Workflow is divided into three steps:
1. Trim adapter and low quality bases from short sequence reads (FASTQ format)
2. Align reads to a reference genome
3. Merge reads originating from the same biological sample
"""

# Config file
configfile: "config.yaml"

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
    
def get_intervals(ref, main_size, backup_size=3000):
    
    """
    Define intervals for scatter-gather variant calling.
    """
    
    intervals = []
    
    with open(ref, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            
            start = 0
            tmp_int = []
            s=record.seq
            scaff_regex = "N"*main_size+ "+"
            intervals += chunk_by_gap(record, main_size, backup_size)
            
#     o = open(config["intervals"], 'w')
#     o.write('\n'.join(intervals))
    return intervals

def chunk_by_gap(record, main_size, backup_size=3000):
    
    """
    Divide reference genome assembly up at gaps of specified size.
    If dividing by gaps of main size fails, dividing at gaps of backup size is attempted
    """
    
    start = 0
    tmp_int = []
    scaff_regex = "N"*main_size
    
    for match in re.finditer(scaff_regex, str(record.seq)):
        tmp_int.append("{}\t{}\t{}".format(record.id, start, start+main_size))
        start = match.end() + 1
    
    # Was it divided up? If not, try a smaller size.
    # Prevents potato chr00 from being run as one chunk, which takes forever
    backup_regex = "N"*backup_size+"+"
    if start == 0:
        for match in re.finditer(backup_regex, str(record.seq)):
            tmp_int.append("{}\t{}\t{}".format(record.id, start, start+backup_size))
            start = match.end() + 1
    else:
        pass
        
    # handle last interval
    tmp_int.append("{}\t{}\t{}".format(record.id, start, len(record)))
    return tmp_int

# Try looking for a scaffold intervals file in data/intervals
# If not available, get_intervals() makes a new intervals file that can be used at restart
try:
    ifh = open(config["intervals"], 'r')
    intervals = []
    for line in ifh:
        intervals.append(line.rstrip())
except FileNotFoundError:
    intervals = get_intervals(config["genome"], config["split_gap_length"])

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
    
# Subfolder wildcard may match: anything ending with a slash, or nothing
# This is to support reference being in top-level or in data/genome/
wildcard_constraints:
    subfolder=".*/|"

# master rule
rule all:
    input:
        ["data/merged/{}.bam".format(x) for x in units.index.levels[0]],
        ["data/summarized_depth/{}_depth_summary.tsv".format(x) for x in units.index.levels[0]],
        config["genome"]+".bwt",
        config["genome"].split(".f")[0] + "_GCN.tsv",
        # config["genome"].split(".f")[0] + "_nongap_{}.bed".format(str(config["split_gap_length"])),
        ["data/medians/{}.bed".format(x) for x in units.index.levels[0]],
        "data/calls/all-calls.vcf"

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
    shell: """
        samtools merge {output} {input}
        samtools index {output}
    """

rule index:
    input:
        "data/merged/{sample}.bam"
    output:
        "data/merged/{sample}.bam.bai"
    shell:
        "samtools index {input}"

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
        win="{{subfolder}}10k_{{ref}}.{ext}.bed".format(ext=config["genome"].rsplit(".")[-1]),
        genome="{{subfolder}}{{ref}}.{ext}".format(ext=config["genome"].rsplit(".")[-1])
    output:
        "{subfolder}{ref}_GCN.tsv"
    shell: """
        python scripts/gcn_content_summarizer.py {input.win} {input.genome} > {output}
    """

# Find all N regions in genome.
# To test: bedtools getfasta -fi ref.fa -bed ref_gaps.bed | less
# rule gap_bed:
#     input:
#         "{{subfolder}}{{ref}}.{ext}".format(ext=config["genome"].rsplit(".")[-1])
#     output:
#         "{subfolder}{ref}_gaps_{gaplength}.bed"
#     shell: """
#         python scripts/gap_finder.py {input} {wildcards.gaplength} > {output}
#     """
# 
# # Apply complement to find all non-N regions in genome.
# # To test: bedtools getfasta -fi ref.fa -bed ref_nongap.bed | grep "N"
# rule safe_bed:
#     input:
#         gaps="{subfolder}{ref}_gaps_{gaplength}.bed",
#         gen="{}.contigs".format(config["genome"])
#     output:
#         "{subfolder}{ref}_nongap_{gaplength}.bed"
#     shell: """
#         bedtools complement -i {input.gaps} -g {input.gen} > {output}
#     """

rule call_variants:
    input:
        ref=config["genome"],
        samples=["data/merged/{}.bam".format(x) for x in units.index.levels[0]],
        indexes=["data/merged/{}.bam.bai".format(x) for x in units.index.levels[0]],
    output:
        "data/calls/intervals/{interval}.vcf"
    params:
        parsed_interval=lambda wildcards: "\t".join(wildcards.interval.rsplit("_", 2)),
        freebayes_options=config["freebayes"]["params"]
    log:
        "log/freebayes/{interval}.log"
    shadow: "full"
    shell: """
        echo {params.parsed_interval} > {wildcards.interval}.bed
        freebayes -t {wildcards.interval}.bed \
            --fasta-reference {input.ref} \
            --bam {input.samples} \
            --vcf {output} \
            {params.freebayes_options} \
            2> {log}
        """

# JVM might need more memory
rule merge_variant_calls:
    input:
        ["data/calls/intervals/{}.vcf".format(re.sub("\t", "_", x)) for x in intervals]
    output:
        "data/calls/all-calls.vcf"
    log:
        "log/merge_vcf.log"
    shell: """
        picard GatherVcfs \
            $(echo {input} | sed 's/data/I=data/g') \
            O={output} \
            2> {log}
        """

rule mosdepth:
    input:
        "data/merged/{sample}.bam"
    output:
        "data/merged/{sample}.mosdepth.global.dist.txt",
        "data/merged/{sample}.mosdepth.summary.txt",
        "data/merged/{sample}.mosdepth.region.dist.txt",
        "data/merged/{sample}.per-base.bed.gz"
    params:
        interval=10000,
        qual=20
    threads: 4
    shell: """
        mosdepth -t {threads} -b {params.interval} -Q {params.qual} {wildcards.sample} {input}
    """
    
################# DRAFT RULES TO SPEED UP DEPTH CALCULATIONS #################
    
rule window_median_depth:
    input:
        "data/merged/{sample}.per-base.bed.gz"
    output:
        "data/merged/{sample}_median_{interval}.bed"
    shell: """
        python scripts/tabix_median.py -w {wildcards.interval} -b {input} -o {output}
    """

rule gather_window_depth:
    input:
        lambda sample: ["data/merged/sample_median_{}.bed".format(re.sub('\t', "_", x) for x in intervals)]
    output:
        "data/medians/{sample}.bed"
    shell: """
        cat {input} > {output}
    """