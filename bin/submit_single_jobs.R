# R-script to guide jobs that can not be run with the pipeline script
#
# An example is for the coverage profile. When running the full pipeline
# for ad-chc project as well as nijmegen, since we do not store the raw data
# the pipeline stops.

# Libraries
library(data.table)
library(rslurm)
library(stringr)

# Input data
# input data is aligned bam file to hg38 (ad-chc project, sequenced in Amsterdam)
aligned_bam = system('ls /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/*.hifi*hg38*bam', intern = T)
# output coverage file where to add lines
outfile = '/project/holstegelab/Share/pacbio/data_processed/coverage_samples.txt'

# Functions
# function to run the coverage profiling
coverage_profile = function(bam_file, prefix, sample_name){
    cmd = paste0("mosdepth -n --fast-mode --by 500 ", prefix, " ", bam_file, " --threads 3; grep -w total ", paste0(prefix, '.mosdepth.summary.txt'), " | sed 's/total/", sample_name, "/g' | cut -f1-4 >> /project/holstegelab/Share/pacbio/data_processed/coverage_samples.txt")
    # add results to main coverage file
    system(paste0("grep -w total ", paste0(prefix, '.mosdepth.summary.txt'), " | sed 's/total/", sample_name, "/g' | cut -f1-4 >> /project/holstegelab/Share/pacbio/data_processed/coverage_samples.txt"))
    res = system(cmd)
    return(res)
}

# Main
# 1. coverage stats -- main loop to find files that need to be processed
    df = data.frame(bam_file = as.character(), prefix = as.character(), sample_name = as.character())
    for (bam in aligned_bam){
        out_file = str_replace_all(bam, '.bam', '.coverage')
        sample_name = str_split_fixed(unlist(strsplit(bam, '/'))[length(unlist(strsplit(bam, '/')))], '.ccs', 2)[, 1]
        # check if files were already done
        if (!file.exists(paste0(out_file, '.mosdepth.summary.txt'))){
            tmp = data.frame(bam_file = bam, prefix = out_file, sample_name = sample_name)
            df = rbind(df, tmp)
        }
    }
    # calculate coverage summary
    sjob <- slurm_apply(coverage_profile, df, jobname = 'test_coverage', cpus_per_node = 3, submit = TRUE)
