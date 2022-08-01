# R-script to guide jobs that can not be run with the pipeline script
#
# An example is for the coverage profile. When running the full pipeline
# for ad-chc project as well as nijmegen, since we do not store the raw data
# the pipeline stops.

# Libraries
library(data.table)
library(rslurm)
library(stringr)

# Functions
# function to run the coverage profiling
coverage_profile = function(bam_file, prefix, sample_name){
    # previous commands using mosdepth
    #cmd = paste0("mosdepth -n --fast-mode --by 500 ", prefix, " ", bam_file, " --threads 3; grep -w total ", paste0(prefix, '.mosdepth.summary.txt'), " | sed 's/total/", sample_name, "/g' | cut -f1-4 >> /project/holstegelab/Share/pacbio/data_processed/coverage_samples.txt")
    #system(paste0("grep -w total ", paste0(prefix, '.mosdepth.summary.txt'), " | sed 's/total/", sample_name, "/g' | cut -f1-4 >> /project/holstegelab/Share/pacbio/data_processed/coverage_samples.txt"))
    # create file where to store coverage information
    cmd = paste0("printf 'GLOBAL_MEDIAN_READ_LENGTH\tGLOBAL_COVERAGE\tMAPPED_COVERAGE\tMAPPED_READS\tALT_COVERAGE\tALT_READS\tUNMAPPED_COVERAGE\tUNMAPPED_READS\n", sample_name, "\t' > ", prefix)
    system(cmd)
    # calculate coverage info
    cmd = paste0("/project/holstegelab/Software/nicco/tools/asbt/build/asbt cov -g 3088000000 ", bam_file, " >> ", prefix)
    system(cmd)
    # add results to main coverage file
    cmd = paste0("tail -1 ", prefix, " >> /project/holstegelab/Share/pacbio/data_processed/coverage_smrt_cells.txt")
    system(cmd)
    return(cmd)
}

# function to copy ccs output to dcache
copy_ccs_dcache = function(ccs_path, smrt_id){
    # define output directory based on project type
    if (length(grep('ad_centenarians', ccs_path)) >0){
        out_dir = 'tape/data_processed/ccs/ad_centenarians/'
    } else if (length(grep('nijmegen', ccs_path)) >0){
        out_dir = 'tape/data_processed/ccs/ad_centenarians/'
    } else if (length(grep('blood', ccs_path)) >0){
        out_dir = 'tape/data_processed/ccs/blood_brain_child/'
    } else {
        out_dir = 'tape/data_processed/ccs/other/'
    }
    
    # set working directory
    tmp = unlist(strsplit(ccs_path, '/'))
    wd = paste(tmp[1:(length(tmp)-1)], collapse = '/')
    setwd(wd)
    
    # then generate md5sum checks file
    outmd5 = str_replace_all(smrt_id, '.ccs.bam', '.md5')
    system(paste0('md5sum ', smrt_id, '* > ', outmd5))

    # copy to dcache
    # ccs.bam
    system(paste0('rclone copy -vv --progress --multi-thread-streams 1 --config /project/holstegelab/Data/dcache.conf ', smrt_id, ' dcache:', out_dir))
    # index
    system(paste0('rclone copy -vv --progress --multi-thread-streams 1 --config /project/holstegelab/Data/dcache.conf ', smrt_id, '.pbi dcache:', out_dir))
    # md5
    system(paste0('rclone copy -vv --progress --multi-thread-streams 1 --config /project/holstegelab/Data/dcache.conf ', outmd5, ' dcache:', out_dir))

    # and check the md5 for a correct copy
    md5check = system(paste0('md5sum -c ', outmd5), intern = T)
    # if the md5 is fine, then we can safely remove the input from spider
    md5check = as.data.frame(str_split_fixed(md5check, ' ', 2))
    if (unique(md5check$V2 == 'OK')){
        # copy was successful, we can safely delete the file
        print('copy was successful, removing file from spider.')
        system(paste0('rm ', ccs_path))
        system(paste0('rm ', ccs_path, '.pbi'))
        system(paste0('rm ', outmd5))
        exitus = 'dcache'
    } else {
        print('!!! copy was NOT successful, keeping the file in spider.')
        exitus = 'spider'
    }
    return(exitus)
}

# Main
# 1. coverage stats -- main loop to find files that need to be processed
    # Input data is aligned bam file to hg38 (ad-chc project, sequenced in Amsterdam)
    aligned_bam = system(paste0('find /project/holstegelab/Share/pacbio/data_processed/nijmegen/ -name "*ccs*.hifi*hg38*bam"'), intern = T)
    # output coverage file where to add lines
    outfile = '/project/holstegelab/Share/pacbio/data_processed/coverage_smrt_cells.txt'
    # create dataframe for submission
    df = data.frame(bam_file = as.character(), prefix = as.character(), sample_name = as.character())
    for (bam in aligned_bam){
        out_file = str_replace_all(bam, '.bam', '.coverage_summary.txt')
        sample_name = str_split_fixed(unlist(strsplit(bam, '/'))[length(unlist(strsplit(bam, '/')))], '.ccs', 2)[, 1]
        # check if files were already done
        if (!file.exists(out_file)){
            tmp = data.frame(bam_file = bam, prefix = out_file, sample_name = sample_name)
            df = rbind(df, tmp)
        }
    }
    # calculate coverage summary
    sjob <- slurm_apply(coverage_profile, df, jobname = 'test_coverage', nodes = 6, cpus_per_node = 8, submit = TRUE)

# 2. copy to dcache ccs output (the biggest file)
    # input data is the very first output of the pipeline (ccs.bam)
    # ad_chc amsterdam
    input_ccs_ad_chc = system('find /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/ -name "m*ccs.bam"', intern = T)
    # blood brain child
    input_ccs_bbc = system('find /project/holstegelab/Share/pacbio/data_processed/blood_brain_child/ -name "m*ccs.bam"', intern = T)
    # nijmegen
    input_ccs_nijmegen = system('find /project/holstegelab/Share/pacbio/data_processed/nijmegen/ -name "m*ccs.bam"', intern = T)
    # anke samples
    input_ccs_anke = system('find /project/holstegelab/Share/pacbio/data_processed/Anke_samples/ -name "m*ccs.bam"', intern = T)
    # combine samples
    all_samples = c(input_ccs_ad_chc, input_ccs_bbc, input_ccs_nijmegen, input_ccs_anke)
    # extract smrt id 
    all_samples_info = data.frame(ccs_path = all_samples, smrt_id = NA, stringsAsFactors = F)
    for (i in 1:nrow(all_samples_info)){
        tmp = unlist(strsplit(all_samples_info$ccs_path[i], '/'))
        all_samples_info$smrt_id[i] = tmp[length(tmp)]
    }
    
    # also collect data that is already in dcache
    already_done = system('find ~/dcache/tape/data_processed/ccs -name "m*ccs.bam"', intern = T)
    # extract smrt id
    already_done_info = data.frame(ccs_path = already_done, smrt_id = NA, stringsAsFactors = F)
    for (i in 1:nrow(already_done_info)){
        tmp = unlist(strsplit(already_done_info$ccs_path[i], '/'))
        already_done_info$smrt_id[i] = tmp[length(tmp)]
    }
    already_done_info$status = 'dcache'

    # exclude the already done ones
    all_samples_info = all_samples_info[which(!(all_samples_info$smrt_id %in% already_done_info$smrt_id)),]

    # add date
    all_samples_info$date = str_split_fixed(all_samples_info$smrt_id, '_', 3)[, 2]
    all_samples_info = all_samples_info[order(all_samples_info$date), ]
    all_samples_info$date = NULL

    # submit job to slurm
    tmp_df = all_samples_info[1:30, ]
    sjob <- slurm_apply(copy_ccs_dcache, tmp_df, jobname = 'copy_ccs_dcache', nodes = 30, cpus_per_node = 1, submit = TRUE)
