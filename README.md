# Oak

Example Snakemake workflow for sequence alignment and variant calling. Handles technical and biological replicates using the [rnaseq-star workflow example](https://github.com/snakemake-workflows/rna-seq-star-deseq2).


## Usage

### STEP 1: Install miniconda and git
If you are using a work or lab server, ask your sysadmin if git and conda are installed already. If so, skip to STEP 2.

To download miniconda:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

The installer will ask you some questions to complete installation. Review and accept the license, accept or change home location, and answer yes to placing it in your path.

To finish configuring miniconda:
```
source $HOME/.bashrc
conda config --add envs_dirs ./.conda/envs
conda config --add pkgs_dirs ./.conda/pkgs
```

To install git:
```
conda install git
```

### STEP 2: Clone the repository

In the terminal, navigate to your preferred location and clone this repository.

```
git clone https://github.com/kramundson/oak.git
cd oak
```

### STEP 3: Build and activate the conda environment
When you __build the conda environment__, Conda obtains all the software listed in `environment.yaml`.
```
conda env create -f environment.yaml
```

Finally, you will need to __activate the environment__. The environment is named "oak," and the software will only be accessible while the environment is active.
```
source activate oak
```

You only need to build the environment once, but if you log out and log in again, you'll need to reactivate the environment before running this workflow.

If you want to deactivate the environment, you can use the command `conda deactivate`.

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

An example units.tsv file is included in the workflow.

### STEP 7: Run Snakemake
To run the workflow with 24 cores and save the printed output to `run_1.out`:
```
snakemake --cores 24 > run_1.out 2>&1
```
