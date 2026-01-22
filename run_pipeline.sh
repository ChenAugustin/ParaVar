# Setup environment
#------------------

CONDA_BASE=$(conda info --base) # Get path to conda
source $CONDA_BASE/etc/profile.d/conda.sh # Import the conda command
conda activate env-snakemake-8.30.0 # Activate Snakemake conda environment

# Setup script and output paths
#------------------------------

snake_dir="/crex/proj/sllstore2017073/private/GRC/Augustin/snakemake/paralog_aware_pipeline/" # Where Snakemake lives
output_dir="/crex/proj/sllstore2017073/private/GRC/Augustin/analysis/snakemake_output-disable_mapping_test/" # Where Snakemake config and output files will be saved

# Create a unique log directory for this run
# Do it now so that we can save logs from the bash commands + a snapshot of the pipeline incl. config etc.
logs_dir="${output_dir}logs-$(date '+%Y%m%d-%H%M%S')/"; mkdir "${logs_dir}" # Create unique log directory, frozen in time
cp -r "${snake_dir}"* "${logs_dir}" # Save snapshot of pipeline
git --git-dir "${snake_dir}.git" rev-parse HEAD > "${logs_dir}git_id" # Save the last commit ID of the git repository
mkdir "${logs_dir}slurm_logs" # Create directory to store --output and --error files from SLURM
mkdir "${logs_dir}rule_logs" # Create directory to store logs from the stderr and stdout of each rule iteration

# Run Snakemake!
#---------------

snakemake --directory "${snake_dir}" \
          --dry-run \
          --rerun-triggers {input,params} \
          --snakefile "${snake_dir}workflow/Snakefile" \
          --printshellcmds \
          --config output_dir="${output_dir}" \
                   logs_dir="${logs_dir}" \
          --profile "${snake_dir}profiles/slurm" \
          --slurm-logdir "${logs_dir}slurm_logs" \
          --slurm-keep-successful-logs \
          1> "${logs_dir}snakemake.out" \
          2> "${logs_dir}snakemake.err"
