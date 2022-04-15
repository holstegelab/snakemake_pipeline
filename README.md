# Snakemake pipeline(s)
Pipeline for the analysis of Pacbio data. The pipeline is subdivided into two main sub-pipelines. The first comprises the analysis of raw data from the sequencer. The second analysis comprises the merging of analyzed data (aligned hifi and non-hifi data) belonging to the same sample.

# Main Pipeline (Snakefile)
## Steps in the pipeline
1. Takes raw sequencing output and run ccs algorithm to generate hifi reads. We retain non-hifi reads, by specifying `--min-passes 0` (number of passes) and `--min-rq 0` (read quality). Kinetics information are also kept (`--all-kinetics`).
2. Run primrose to generate methylation profiles for each read
3. Split hifi reads from non-hifi reads
4. Alignment of hifi and non-hifi reads to GRCh38 (hg38)
5. Alignment of hifi and non-hifi reads to CHM13
6. Perform sample check by comparing pacbio genotypes with GWAS array genotypes

## To run the pipeline, use the following commands:
### It is adviseable to run the pipeline in a screen process
`screen -S name_screen_process`

### Load conda environment py37 to use snakemake
`conda activate py37`

### Main command
`snakemake -s path/to/Snakefile --latency-wait 60 --configfile path/to/config_XXX.yml --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config path/to/config_cluster.yml --jobs 1`

### Parameters
The only parameter to define is the configuration file (`--configfile`). This is a plain-text file containing the directory of input files and the desider directory of output files. An example config file should look like:

`IN_DIR : "path/to/input/directory"`

`OUT_DIR : "path/to/desider/output/directory"`

# Second Pipeline (Snakefile_merge_and_pileup)
## Steps in the pipeline
1. Merge and index hifi data aligned to GRCh38
2. Merge and index non-hifi data aligned to GRCh38
3. Merge and index hifi data aligned to chm13
4. Merge and index non-hifi data aligned to chm13
5. Run pileup analysis for methylation profiles for GRCh38-aligned (merged) data
6. Run pileup analysis for methylation profiles for chm13-aligned (merged) data
7. Run deepvariant analysis on GRCh38-aligned data

## To run the pipeline, use the following commands:
### It is adviseable to run the pipeline in a screen process
`screen -S name_screen_process`

### Load conda environment py37 to use snakemake
`conda activate py37`

### Main command
`snakemake -s path/to/Snakefile_merge_and_pileup --latency-wait 60 --configfile path/to/config_merge/config_XXX.yml --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config path/to/config_cluster_merge_and_pileup.yml --jobs 1`

### Parameters
Similarly to the previous pipeline, also this takes a configuration file as input. This should be a plain-text file containing the directory and file prefix of input and output files. An example config file could look like this:

`IN_FILES : "path/to/c1_blood/m64037e_210709_135828,path/to/c1_blood/m64050_201127_132810,path/to/c1_blood/m64050_210223_091925,path/to/c1_blood/m64050_210709_135541"`

`OUT_FILE : "path/to/c1_blood/c1_merged"`