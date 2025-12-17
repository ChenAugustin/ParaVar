#### Install plugin
To execute the paralog_aware workflow on a SLURM cluster, you will need to install the [SLURM executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) into the paralog_aware environment.
```shell
conda activate paralog_aware
mamba install -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```
