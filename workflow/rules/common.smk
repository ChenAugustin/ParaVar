##############################
# Import Python packages
##############################

import sys
import os
import pandas as pd


##############################
# Import variables
##############################

# Store config values in variables for:
#     - Snakemake routines (only to enhance clarity)
#     - Shell within Snakemake rules (mandatory)

# Essentials
#===========

# Output paths
OUTPUT_DIR = config["output_dir"]
LOGS_DIR = config["logs_dir"]

# Standard directory structure for rule output
RES_DIR_RULE = OUTPUT_DIR + "results/"
LOG_DIR_RULE = LOGS_DIR + "rule_logs/"

# Reference genome
REFGENOME = config["refGenome"]
REFPATH = config["refPath"]

# Settings
#=========

# Trimmomatic
READ_TYPE = config["read_type"]
TRIMMOMATIC_SETTING = config["trimmomatic_setting"]

# Site calling
#-------------

# Intervals to parallelize
CONTROLSEQ = config["controlseq"]
N_INTERVALS_CONTROLSEQ = config["n_intervals_controlseq"]
N_DIGITS_CONTROLSEQ = len(str(N_INTERVALS_CONTROLSEQ))
INTERVALS_CONTROLSEQ = [f"{x:0{N_DIGITS_CONTROLSEQ}d}" for x in range(N_INTERVALS_CONTROLSEQ)]

TARGETSEQ = config["targetseq"]
N_INTERVALS_TARGETSEQ = config["n_intervals_targetseq"]
N_DIGITS_TARGETSEQ = len(str(N_INTERVALS_TARGETSEQ))
INTERVALS_TARGETSEQ = [f"{x:0{N_DIGITS_TARGETSEQ}d}" for x in range(N_INTERVALS_TARGETSEQ)]

# Site filtering
#---------------

# Basic filters
MIN_MISSING = config["min_missing"]
MIN_QUAL = config["min_qual"]
MIN_MQ = config["min_mq"]
MAX_MQ0F = config["max_mq0f"]

# Mismapping controllib targetseq
SHUF_N_LINES_PER_INTERVAL = config["shuf_n_lines_per_interval"]
SAMPLE_SHEET_PATH = config["samples"]
FRACTION_TARGETSEQ_IN_TARGETLIB = config["fraction_targetseq_in_targetlib"]
VCF_HEADER_MISMAP = config["vcf_header_mismap"]
MAX_MMMCT = config["max_mmmct"]

# Depth filter
TARGETSEQ_SINGLECOPY = config["targetseq_singlecopy"]
MIN_GDEPTH = config["min_gdepth"]
FACTOR_MAX_GDEPTH = config["factor_max_gdepth"]
VCF_HEADER_GDEPTH = config["vcf_header_gdepth"]
MIN_NDP = config["min_ndp"]

# Nucleotide diversity
PIXY_POPULATION_FILE = config["pixy_population_file"]
PIXY_WINDOW_SIZE = config["pixy_window_size"]

# Plot
CIRCOS_KARYOTYPE_TARGETSEQ = config["circos_karyotype_targetseq"]
CIRCOS_CONF_TARGETLIB_TARGETSEQ = config["circos_conf_targetlib_targetseq"]


##############################
# Import essential data
##############################

# Dataframe of samples
def parse_sample_sheet(sample_sheet_path: str) -> pd.DataFrame:
    samples = (
        pd.read_table(sample_sheet_path, sep=",", dtype=str)
        .replace(" ", "_", regex=True)
        .infer_objects(copy=False) # needed to maintain same behavior in future pandas versions
    )
    return samples
samples = parse_sample_sheet(config["samples"])

# Add new columns with the expected paths to trimmed FASTQ files
samples["fq1.trimmomatic_1P"] = samples.apply(lambda row: f"{RES_DIR_RULE}trimmomatic/{row.fastq_id}_1P.fastq.gz", axis=1)
samples["fq2.trimmomatic_2P"] = samples.apply(lambda row: f"{RES_DIR_RULE}trimmomatic/{row.fastq_id}_2P.fastq.gz", axis=1)

# Convert specific columns into lists over which rules can iterate
fastq_ids = samples["fastq_id"].tolist()
sample_ids_controllib = samples.loc[samples["targetlib"] == "FALSE", "sample_id"].unique().tolist()
sample_ids_targetlib = samples.loc[samples["targetlib"] == "TRUE", "sample_id"].unique().tolist()


##############################
# Define functions for rules
##############################

# Generate list of files requested by rule all
def get_output():

    # Initialize empty list
    out = []

    # Request mandatory output
    #-------------------------

    out.append(RES_DIR_RULE + "gather_vcfs_filtered_targetlib_targetseq/targetseq.nindel_filtered.vcf.gz")

    # Request optional output
    #------------------------

    # FastQC/MultiQC output
    if config["run_fastqc"]:
        out.extend(expand(RES_DIR_RULE + "fastqc/{fastq_id}_R1_fastqc.html", fastq_id = fastq_ids))
        out.extend(expand(RES_DIR_RULE + "fastqc/{fastq_id}_R1_fastqc.zip", fastq_id = fastq_ids))
        out.extend(expand(RES_DIR_RULE + "fastqc/{fastq_id}_R2_fastqc.html", fastq_id = fastq_ids))
        out.extend(expand(RES_DIR_RULE + "fastqc/{fastq_id}_R2_fastqc.zip", fastq_id = fastq_ids))
        out.append(RES_DIR_RULE + "multiqc/multiqc_report.html")
    if config["run_fastqc_post"]:
        out.extend(expand(RES_DIR_RULE + "fastqc_post/{fastq_id}_1P_fastqc.html", fastq_id = fastq_ids))
        out.extend(expand(RES_DIR_RULE + "fastqc_post/{fastq_id}_1P_fastqc.zip", fastq_id = fastq_ids))
        out.extend(expand(RES_DIR_RULE + "fastqc_post/{fastq_id}_2P_fastqc.html", fastq_id = fastq_ids))
        out.extend(expand(RES_DIR_RULE + "fastqc_post/{fastq_id}_2P_fastqc.zip", fastq_id = fastq_ids))
        out.append(RES_DIR_RULE + "multiqc_post/multiqc_report.html")

    # Statistics
    out.append(RES_DIR_RULE + "pixy_nucleotide_diversity_filtered_targetlib_targetseq/targetseq.nindel_filtered_pi.txt")

    # Plots
    # out.append(RES_DIR_RULE + "circos_plot_stats_filtered_targetlib_targetseq/circos.svg")
    out.append(RES_DIR_RULE + "rideogram_plot_stats_filtered_targetlib_targetseq/scorability.svg")
    out.append(RES_DIR_RULE + "rideogram_plot_stats_filtered_targetlib_targetseq/pi.svg")

    return out


# If input is short paired-end reads, parses samples pd.DataFrame and extract fq1 and fq2 given wc.fastq_id for Trimmomatics
def get_reads_to_trim(wc):

    # Get row that contains wc.fastq_id
    row = samples.loc[samples["fastq_id"] == wc.fastq_id]

    # Extract fq1 and fq2 from the row
    fq1 = row["fq1"].values[0]
    fq2 = row["fq2"].values[0]

    # Return a dictionnary that can be directly used as Snakemake `input:`
    if config["read_type"] == "short_pe":
        return {"fq1": fq1, "fq2": fq2}
    elif config["read_type"] == "long_hifi":
        return {"fq1": fq1}


# Get correct set of input reads depending on config["run_trimmomatic"] and config["read_type"]
def get_reads_to_map(wc):

    # Get row(s) that contains wc.sample_id
    df = samples[samples.sample_id == wc.sample_id]

    # Extract list of FASTQ from wc.sample_id to map
    #
    # If use HiFi long reads as input, then return a single FASTQ file
    # Note that HiFi long reads don't need trimming
    if config["read_type"] == "long_hifi":
        return {"r1": df["fq1"].tolist()}
    # Else, return two FASTQ files as paired-end short reads by default
    else:
        # If Trimmomatic is run, then get trimmed reads from columns fq1.trimmomatic_1P/fq2.trimmomatic_2P
        # Else, get reads from columns fq1/fq2
        if config["run_trimmomatic"]:
            return {"r1": df["fq1.trimmomatic_1P"].tolist(),
                    "r2": df["fq2.trimmomatic_2P"].tolist()}
        else:
            return {"r1": df["fq1"].tolist(),
                    "r2": df["fq2"].tolist()}
