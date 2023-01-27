# Snakemake pipeline(s)
Pipeline for the analysis of Pacbio data. The pipeline is subdivided into two main sub-pipelines. The first comprises the analysis of raw data from the sequencer. The second analysis comprises the merging of analyzed data (aligned hifi and non-hifi data) belonging to the same sample.

# Main Pipeline (Snakefile)
The latest development of the pipeline file is able to automatically (1) retrieve the data that need to be processed, (2) write the necessary configuration files, and (3) submit the pipeline.
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
### It is adviseable to run the script interactively, to avoid the submission of jobs using data that is still actively being copied to dcache. 
The script to use is `run_snakemake_pipeline.py`. This script will open a new screen window with a defined name, and will run the pipeline within that screen window. Alternatively (yet not recommended), if you are sure that all data in dcache has been copied correctly, you can also do:
`python3 run_snakemake_pipeline.py`
As a consequence, all files present in dcache that have not been analyzed yet, will be analyzed.

### Main command
Behind the scenes, the command that is actually use to run the pipeline is contained in `submit_copy_and_pipeline.sh`. This bash script takes a single argument as parameter, that is, the configuration file containing the information of the sample to be processed. This bash script calls the actual `snakemake` command that guides the pipeline.

### Parameters
The only parameter to define to the bash script is the configuration file containing information of the sample (smrt cell) to be analyzed. Briefly, a configuration file contains the path to the data to be copied from dcache, and the output directory of the analyzed data. An example config file should look like:
`IN_DIR : "tape/path/to/subreads/mXXXXX_XXXXXX_XXXXX.subreads.bam"`

`OUT_DIR : "path/to/desider/output/directory"`
The `IN_DIR` parameter links to dcache storage. The path to dcache storage is `~/dcache`.

### Updates
Due to the recent update to SMRTlink v11, the previous CCS algorithm version 6.0 did not work anymore. As a consequence, an updated version of the CCS algorithm (version 6.4) was needed, and it is now installed in conda environment smrtlink11. Because only the machines in Nijmegen have (for now) been updated to SMRTlink v11, at the moment both the previous pipeline (based on CCS v6) and the new pipeline (based on CCS v6.4) are present. The updated scripts are `submit_copy_and_pipeline_updated.sh` and `Snakefile_with_copy_update.py`.

# Second Pipeline (Snakefile_merge_and_pileup)
Similar to the first pipeline, the latest development of this pipeline is also able to automatically (1) retrieve the data that need to be processed, (2) write the necessary configuration files, and (3) submit the pipeline.
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
### Althought it is adviseable to run the script interactively, the script can be run directly. 
The script to use is `run_snakemake_merge.py`. This script will open a new screen window with a defined name, and will run the pipeline within that screen window. `python3 run_snakemake_merge.py`
At the moment, this script will merge only samples for which the combined coverage is >12x.

### Main command
Behind the scenes, the command that is actually use to run the pipeline is contained in `submit_merge_and_pileup.sh`. This bash script takes a single argument as parameter, that is, the configuration file containing the information of the sample to be processed. This bash script calls the actual `snakemake` command that guides the pipeline.

### Parameters
The only parameter to define to the bash script is the configuration file containing information of the sample (smrt cell) to be analyzed. Briefly, a configuration file contains the path to the datasets to be combined together, and the output directory of the analyzed data. An example config file should look like:
`IN_FILES : "path/to/sample_1/m64037e_210709_135828,path/to/sample_2/m64050_201127_132810,path/to/sample_3/m64050_210223_091925,path/to/sample_4/m64050_210709_135541"`
`OUT_FILE : "path/to/sample_merged/sample_merged"`
In the command above, 4 samples will be combined together. The directory and prefix of the combined data is defined with `OUT_FILE`.

# Additional scripts
## submit_single_jobs.R
Script in R that run a single operation on a set of .bam files. At the moment, the only operation is the coverage calculation. This was necessary for the files for which input raw data is deleted after the pipeline. Now the coverage operation is implemented as part of the pipeline.
A second part of this script is still used to copy raw CCS output to dcache as this files are the largest (in terms of storage).

## generate_freeze_coverage.R
Script in R that takes the summary statistics of the coverage of each sample, merge them together and then add information such as the project (blood-brain-child, ad-centenarians, anke), sequencing center (vumc, nijmegen) as well as sample information (sample id, phenotype, data path).
This script should be run every week to generate a new freeze. The results should then be put in the file on the researchdrive at `pacbio/YYYYMMDD_pacbio_sequencing_Nicco.xlsx`.

