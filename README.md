# Snakemake pipeline(s)
Pipeline for the analysis of Pacbio data. The pipeline is subdivided into two main sub-pipelines. The first comprises the analysis of raw data from the sequencer. The second analysis comprises the merging of analyzed data (aligned hifi and non-hifi data) belonging to the same sample.

# Main Pipeline (Snakefile)
## Steps in the pipeline
1. Copy of sequencing data from dcache to local storage
2. MD5 checksum check
3. CCS generation with custom parameters: `--min-passes 0` (number of passes), `--min-rq 0` (read quality), and `--all-kinetics` (kinetics information are kept)
2. Run primrose to generate methylation profiles for each read
3. Split hifi reads from non-hifi reads
4. Alignment of hifi and non-hifi reads to GRCh38 (hg38)
5. Alignment of hifi and non-hifi reads to CHM13
6. Perform sample check by comparing pacbio genotypes with GWAS array genotypes
7. Coverage summary per chromosome and total

## To run the pipeline, use the following commands:
### It is adviseable to run the pipeline in a screen process
`screen -S name_screen_process`

### Load conda environment py37 to use snakemake
`conda activate py37`

### Main command
`snakemake -s path/to/Snakefile_with_copy.py --latency-wait 60 --configfile path/to/config_XXX.yml --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config path/to/config_cluster.yml --jobs 1`

### Parameters
The only parameter to define is the configuration file (`--configfile`). This is a plain-text file containing the directory of input files and the desider directory of output files. An example config file should look like:

`IN_DIR : "tape/path/to/subreads/mXXXXX_XXXXXX_XXXXX.subreads.bam"`

`OUT_DIR : "path/to/desider/output/directory"`

The `IN_DIR` parameter links to dcache storage. The path to dcache storage is `~/dcache`.

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
`conda activate cpg`

### Main command
`snakemake -s path/to/Snakefile_merge_and_pileup --latency-wait 60 --configfile path/to/config_merge/config_XXX.yml --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask}" --cluster-config path/to/config_cluster_merge_and_pileup.yml --jobs 1`

### Parameters
Similarly to the previous pipeline, also this takes a configuration file as input. This should be a plain-text file containing the directory and file prefix of input and output files. An example config file could look like this:

`IN_FILES : "path/to/sample_1/m64037e_210709_135828,path/to/sample_2/m64050_201127_132810,path/to/sample_3/m64050_210223_091925,path/to/sample_4/m64050_210709_135541"`

`OUT_FILE : "path/to/sample_merged/sample_merged"`

In the command above, 4 samples will be combined together. The directory and prefix of the combined data is defined with `OUT_FILE`.

# Additional scripts
## submit_single_jobs.R
Script in R that run a single operation on a set of .bam files. At the moment, the only operation is the coverage calculation. This was necessary for the files for which input raw data is deleted after the pipeline. Now the coverage operation is implemented as part of the pipeline.

## generate_freeze_coverage.R
Script in R that takes the summary statistics of the coverage of each sample, merge them together and then add information such as the project (blood-brain-child, ad-centenarians, anke), sequencing center (vumc, nijmegen) as well as sample information (sample id, phenotype, data path).
This script should be run every week to generate a new freeze. The results should then be put in the file on the researchdrive at `pacbio/YYYYMMDD_pacbio_sequencing_Nicco.xlsx`.

## submit_merge_and_pileup.sh
Bash script that is called by another script (`run_snakemake_merge.py`) in order to submit the snakemake pipeline to merge single smrt-cells results together, do the pileup analysis and deepvariant.

## run_snakemake_merge.py
This script can be used to merge data processed of single smrt cells together. At the moment, only samples with a combined coverage >15 are merged. At this moment (19 July 2022), 70 samples satisfy these parameters and are being merged. This script performs:
1. collection of all single smrt cells results
2. collection of their likely samples
3. calculation of total coverage after combining smrt cells from the same sample
4. if the combined total coverage is larger than 15, a configuration file readable by snakemake is created
5. Finally, the `Snakefile_merge_and_pileup.py` pipeline is used to merge results. The submission automatically opens a screen process named `merge_{SAMPLE_NAME}`, loads the required conda environment `conda activate cpg` and finally runs the `submit_merge_and_pileup.py` script. This latter takes only 1 argument, that is, the configuration file needed by snakemake.
`run_snakemake_merge.py` generates (when run for the first time) and then updates 2 files:
1. a file containing the config_files that have been submitted through the pipeline (at `/project/holstegelab/Software/snakemake_pipeline/config/config_merge/merged_submitted.txt`)
2. a summary file containing the date of submission, sample identifier and diagnosis, number of smrt cells involved, directory of outputs, and combined coverage (at `/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt`). This latter file will then be copied to the main Excel file. Idea is to run all merging operations on Friday, update the freeze and finally update the Excel file.