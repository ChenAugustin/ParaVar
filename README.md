# ParaVar

![image-flowchart](https://github.com/ChenAugustin/ParaVar/blob/main/docs/img/flowchart.png)

**ParaVar** retrieves SNPs from **tissue-specific genomic regions** (e.g., sex chromosomes, B chromosomes, germline-restricted chromosomes), which are often challenging to analyze due to the presence of numerous paralogs, sometimes sharing near-identical sequences. It does so by mapping a library that contains the tissue-specific "target" regions and one that doesn't, then comparing the two using diverse statistics.

ParaVar is delivered as a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html#) pipeline for enhanced reproducibility and scalable analyses. It produces a standard VCF file as output, ensuring compatibility with most downstream analyses.

<br/>

## Citing ParaVar

If you use ParaVar, please cite:
- [...]
- Also, make sure to cite the tools you used within ParaVar.

<br/>

## Table of contents

* [Data required](#data-required)
* [Installation](#installation)
* [Configs](#configs)
    * [1. Sample sheet](#1.-sample-sheet)
    * [2. Pixy file](#2.-pixy-file)
    * [3. Config file](#3.-config-file)
    * [4. Run file](#4.-run-file)
    * [5. Cluster file](#5.-cluster-file)
* [Run ParaVar](#run-paravar)
 
<br/>

## Data required

- **Target libraries** containing reads from the target regions
- Paired with each target library, a **control library** without reads from the target regions
- A **sample sheet** with pairing information of the target and control libraries
- A reference genome with **coordinates** of the control regions and the target regions
- Coordinates of a sequence in the target regions believed to be in **single copy**

<br/>

## Installation

All you need to start are two dependencies: **conda** and **Snakemake**.

The goal is to install Snakemake in its own separate environment using conda. For this:
1. Install conda via the [Miniforge distribution](https://github.com/conda-forge/miniforge). Follow the instructions in the link.

2. Once conda is installed, use it to create a separate environment in which you install Snakemake. Note that mamba is just a faster version of conda and both are included when installing Miniforge.
   ```bash
   mamba create -c conda-forge -c bioconda -n env-snakemake-8.30.0 snakemake=8.30.0
   ```
3. Finally, install the plugin that handles the job scheduler of your cluster. Here's an example for SLURM:
   ```bash
   mamba install -c conda-forge -c bioconda -n env-snakemake-8.30.0 snakemake-executor-plugin-slurm=1.0.1
   ```
   For other plugins, refer to [Snakemake's documentation](https://snakemake.github.io/snakemake-plugin-catalog/)

Now, all that's left to do is clone the ParaVar GitHub repo:
```bash
git clone https://github.com/ChenAugustin/ParaVar.git
cd ParaVar
```

<br/>

### Version
ParaVar was successfully tested using:
- Miniforge version 24.11.3-2
- Snakemake version 8.30.0
- snakemake-executor-plugin-slurm version 1.0.1

<br/>

## Configs

There are **5 configuration files**, all mandatory.

<br/>

### 1. Sample sheet

Located in `config/samples.csv`, it indicates where to find the data and how it is organized. It must be comma-separated and contain the following columns:

|Field|Description|Examples|
|-----|-----------|-------|
|`targetlib`|Indicates if the library contains the target region|`TRUE`, `FALSE`|
|`individual_id`|Must be the same for the target and control libraries from the same pair|`id1`, `id2`|
|`sample_id`|Must be different for the target and control libraries from the same pair|`id1-control`, `id1-target`|
|`fastq_id`|Must be different for each pair of FASTQ files and correspond to the string before `_R1` and `_R2`|`id1-control-fastq1`|
|`fq1`|Path to the 1st FASTQ file. Name format is strict (see example). Symbolic links with the correct format are supported|`/path/to/id1-control-fastq1_R1.fastq.gz`|
|`fq2`|Path to the 2nd FASTQ file. Name format is strict (see example). Symbolic links with the correct format are supported|`/path/to/id1-control-fastq1_R2.fastq.gz`|
|`bam`|Optional: required if mapping is **disabled** (`config["run_mapping"]` is `false`). Path to the user-provided BAM file|`/path/to/sample_id.bam`|
|`library_size`|Optional: required if estimating the fraction of target reads in target libraries by calculating median read depth across the control region is **disabled** (`config["run_samtools_depth_lib_controlseq"` is `false`). In this case, the fraction of target reads is assumed to be the same across all target libraries. This is required to normalize the fraction of control reads from a target library mismapped to each target site (field `FORMAT/MMCT` in VCF file). Any statistic reflecting overall library size can be used (e.g., number of reads)|`776818512`|

Each row corresponds to a pair of FASTQ files. If several FASTQ files map to the same sample (e.g., multiplexed, multiple sequencing runs, multiple libraries), each FASTQ pair should be provided in a separate row, with the **same `sample_id`**, but **different `fastq_id`**.

<br/>

### 2. Pixy file

Located in `config/pixy_population_file.txt`, it is a tiny file required by **pixy**, the tool used inside ParaVar to compute nucleotide diversity ($\pi$). Contains only 2 columns:

|Field|Description|Examples|
|-----|-----------|-------|
|`sample_id`|Must be the same values as in `config/samples.csv`. Provide only the values of the target libraries as ParaVar's focus is to compute $\pi$ of target libraries across scorable sites of the target region|`id1-target`, `id2-target`|
|`population_id`|Which population each target sample comes from. Required by pixy. If all target samples are from the same population, provide the same string for all|`pop1`|

<br/>

### 3. Config file

Located in `config/config.yaml`, it contains ParaVar's main variables: which step to run, how to filter sites, tool settings, etc.
A detailed description of each variable is provided in `config/config.yaml`

<br/>

### 4. Run file

Located in `run_pipeline.sh`, it is the executable `bash` script required to start ParaVar (see how to use it in **Run ParaVar**). It contains the following key variables:

|Field|Description|Examples|
|-----|-----------|-------|
|`snake_dir`|Directory where ParaVar is located|`path/to/ParaVar`|
|`output_dir`|Directory where output will be saved. Note that each time you run `run_pipeline.sh` with the same `output_dir`, Snakemake will attempt to resume the analysis (if it broke at some step), but a new timestamped log file will be created|`path/to/output_dir`|

The run file also contains the actual `snakemake` command, which you can edit using Snakemake's arguments, such as `--dry-run`

<br/>

### 5. Cluster file

Located in `profiles/slurm/config.yaml`, it contains information on the compute resources to request for each step of the pipeline. A detailed description of each variable is provided in `profiles/slurm/config.yaml`

<br/>

## Run ParaVar

Now that everything is set, let's run ParaVar! As it may take a while for all jobs to complete, we recommend running ParaVar within a `tmux` terminal:
```
tmux new-session -s ParaVar
```
Run ParaVar within the `ParaVar` tmux terminal: 
```
bash run_pipeline.sh
```
Press `Ctrl+b` then `d` to exit the tmux terminal without terminating it

We highly recommend first running `snakemake` after adding `--dry-run` inside `run_pipeline.sh`, as it shows what Snakemake intends to do in `${logs_dir}/snakemake.out`. This allows early troubleshooting and avoids re-running steps when resuming an analysis
