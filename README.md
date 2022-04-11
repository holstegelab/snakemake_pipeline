# snakemake_pipeline
Pipeline for the analysis of Pacbio data

# Steps in the pipeline
1. Takes raw sequencing output and run ccs algorithm to generate hifi reads. We retain non-hifi reads, by specifying --min-passes 0 (number of passes) and --min-rq 0 (read quality). Kinetics information are also kept (--all-kinetics).
2. Run primrose to generate methylation profiles for each read
3. Split hifi reads from non-hifi reads
4. Alignment of hifi and non-hifi reads to GRCh38 (hg38)
5. Alignment of hifi and non-hifi reads to CHM13
6. Perform sample check by comparing pacbio genotypes with GWAS array genotypes

# To run the pipeline, use the following commands:
## It is adviseable to run the pipelin in a screen process
screen -S name_screen_process

## Load conda environment py37 to use snakemake
conda activate py37

## Main command
snakemake -s /project/holstegelab/Software/snakemake_pipeline/bin/Snakefile --latency-wait 60 --configfile /project/holstegelab/Software/snakemake_pipeline/bin/config/config_files_snakemake/config_XXX.yml --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config /project/holstegelab/Software/snakemake_pipeline/config/config_cluster.yml --jobs 1

## Parameters
The only parameter to define is the configuration file (--configfile). This is a plain-text file containing the directory of input files and the desider directory of output files. An example config file should look like:

IN_DIR : "path/to/input/directory"

OUT_DIR : "path/to/desider/output/directory"