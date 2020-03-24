# Staples

FASTQ to BAM to VCF. A bioinformatics staple.

This workflow for sequence alignment and variant calling is implemented with [Snakemake](https://snakemake.readthedocs.io/en/stable/) for ease-of-use and reproducibility. Handles technical and biological replicates following the [rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2) workflow example.

If this is your first time using a Snakemake workflow, you may find it helpful to run our [tutorial](https://github.com/lisakmalins/Snakespeare) first.


## Usage

### STEP 1: Install miniconda and git
If you are using a work or lab server, ask your sysadmin if git and conda are installed already. If so, skip to STEP 2.

To __download miniconda__:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

The installer will ask you some questions to complete installation. Review and accept the license, accept or change home location, and answer yes to placing it in your path.

To finish configuring miniconda:
```
source $HOME/.bashrc
```

To __install git__:
```
conda install git
```

### STEP 2: Clone the repository

In the terminal, navigate to your preferred location and __clone this repository__.

```
git clone https://github.com/kramundson/staples.git
cd staples
```

### STEP 3: Build and activate the conda environment
When you __build the conda environment__, Conda obtains all the software listed in `environment.yaml`.
```
# Recommended: prevent conda from crashing if home folder is not writable
conda config --add envs_dirs ./.conda/envs
conda config --add pkgs_dirs ./.conda/pkgs

# Build staples conda environment
conda env create -f environment.yaml
```

Finally, you will need to __activate the environment__.
```
source activate staples
```

You only need to build the environment once. However, you'll need to activate the environment each time you log in. To deactivate the environment, use the command `conda deactivate`.

### STEP 4: Download the reference genome and reads
You'll need to place the reference genome assembly for your organism in `data/genome`.
- If you are downloading the assembly from the Internet, you can use `wget` to download it directly to the folder.

You will also need to place the reads into the repository.
- If you have `fastq` or `fastq.gz` files already, you may put them directly in `data/reads`.
- If you are downloading from NCBI SRA, you may place the sra files in `data/sra` and the workflow will convert them to fastq.

### STEP 5: Customize `config.yaml`
The config file `config.yaml` allows you to customize the workflow to your needs. An example is included; you will need to edit the config file to run your own data.

### STEP 6: Customize `units.tsv`
The purpose of `units.tsv` is to associate sample information with each of your fastq files. It is a tab-delimited file with four required columns:
- `sample`: biological sample
- `unit`: unique combination of biolgical sample, library and sequencing run
- `fq1`: forward reads fastq file
- `fq2`: reverse reads fastq file

An example `units.tsv` file is included in the workflow.

### STEP 7: Run Snakemake
To run the workflow with 24 cores and save the printed output to `run_1.out`:
```
snakemake --cores 24 > run_1.out 2>&1
```
