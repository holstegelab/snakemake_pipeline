# Snakemake pipeline(s)
Pipeline for the analysis of Pacbio data. The pipeline is subdivided into two main sub-pipelines:
1. The first (1_SMRT_level) comprises the analysis of a single SMRT cell after sequencing
2. The second (2_Sample_level) involves merging data from multiple SMRT cells of the same individual

# SMRT-level
The pipeline was originally develop to automatically download data from tape and process data. Small changes to Snakemake file (line 136) should be made to modify the script accordingly.
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
### Start a new screen window
`screen -S snakemake_job_1`
### Activate conda environment
`conda activate py37`
### Run snakemake
snakemake -s [path/to/snakemake_smrt.py] --latency-wait 60 --configfile [path/to/config_example.yml] --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask} --time {cluster.time}" --cluster-config [path/to/config_cluster.yml] --jobs 1

### Run configuration file
There is only 1 configuration file in which the user needs to define the input subread.bam and the output directory. An example config file is provided and should look like:
`IN_DIR : "tape/path/to/subreads/mXXXXX_XXXXXX_XXXXX.subreads.bam"`
`OUT_DIR : "path/to/desider/output/directory"`
Normally, `IN_DIR` parameter links to the path in the tape where subreads.bam data is stored. It should be sufficient to change this to the data location (along with small tuning in the main snakemake file) to make the pipeline to work.

### Cluster configuration file
The snakemake pipeline is implemented to work with SLURM-based systems. This means that wach job in the pipeline is submitted to a computing cluster using SLURM. A configuration file for this is needed, in order to provide the resources for each submitted job. An example config file is provided. Assuming a job named `ccs`, one can define resources to use with:
`ccs:
  ntasks: 1
  ncpupertask: 50
  time: 2880
`
which defines 1 task, 50 CPUs and 2880 minutes (48 hours)

# Sample-level
This pipeline, given a list of SMRT cells to be combined, performs merging, de-novo assembly, variant calling, copy of the raw data to tape, and check of these copied data.
## Steps in the pipeline
1. Merge and index hifi data aligned to GRCh38
2. Merge and index non-hifi data aligned to GRCh38
3. Merge and index hifi data aligned to chm13
4. Merge and index non-hifi data aligned to chm13
5. Run pileup analysis for methylation profiles for GRCh38-aligned (merged) data
6. Run pileup analysis for methylation profiles for chm13-aligned (merged) data
7. Run deepvariant analysis on GRCh38-aligned data
8. Convert merged hifi to fasta format
9. De-novo assembly using hifi reads and hifiasm (phased output)
10. De-novo assembly using hifi reads and hifiasm (phases are collapsed)
11. De-novo assembly using hifi reads and flye
12. Alignment of de-novo assemblies (hifiasm and flye) to GRCh38
13. Alignment of de-novo assemblies (hifiasm and flye) to chm13
14. BUSCO analysis for genome assembly completeness (based on the different assemblies)
15. Structural variant analysis using pbsv
16. Structural variant analysis using sniffles
17. Copy of the single-run data to dcache
18. Check whether data is correctly copied and move data to trashbin

## To run the pipeline, use the following commands:
### Start a new screen window
`screen -S snakemake_job_merge_1`
### Activate conda environment
`conda activate cpg`
### Run snakemake
snakemake -s [path/to/snakemake_sample.py] --latency-wait 60 --configfile [path/to/config_example_sample.yml] --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask} --time {cluster.time}" --cluster-config [path/to/config_cluster_sample.yml] --jobs 1

### Run configuration file
There is only 1 configuration file in which the user needs to define the input subread.bam(s) and the output directory. An example config file is provided and should look like:
`IN_DIR : "tape/path/to/subreads/smrt_prefix_1,tape/path/to/subreads/smrt_prefix_2"`
`OUT_DIR : "path/to/desider/output/directory/sample_prefix"`

### Cluster configuration file
Similarly to the first step of the pipeline, a cluster configuration file should be defined with the resources to use for each job.