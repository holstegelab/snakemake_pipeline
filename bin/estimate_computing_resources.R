# SCRIPT TO ESTIMATE SINGLE DATA-PROCESSING RESOURCES
# DESCRIPTION
# AT THE MOMENT, THE PIPELINE IS MADE OF THESE STEPS:
# 1. COPY FROM DCACHE
# 2. CCS
# 3. PRIMROSE
# 4. EXTRACT-CCS
# 5. ALIGN CCS - HG38
# 6. ALIGN CCS - CHM13
# 7. ALIGN NON-CCS - HG38
# 8. ALIGN NON-CCS - CHM13
# 9. SAMPLE CHECK
# 10. COVERAGE STAT
# 11. MERGE CCS - HG38
# 12. MERGE NON-CCS - HG38
# 13. MERGE CCS - CHM13
# 14. MERGE NON-CCS - CHM13
# 15. PILEUP - HG38
# 16. PILEUP - CHM13
# 17. DEEPVARIANT

# LIBRARY
    library(ggplot2)
    library(ggbreak)
    library(grid)
    library(vioplot)
    library(plyr)
    library(egg)
    library(ggpubr)

# FUNCTIONS
    # function to translate time (m, d, h) into hour time
    transl = function(tmp_res){
        tmp_res = stringr::str_replace_all(tmp_res, 'ms', '/3600/1000+')
        tmp_res = stringr::str_replace_all(tmp_res, 'us', '/3600/1000000+')
        tmp_res = stringr::str_replace_all(tmp_res, 'ns', '/3600/1000000000+')
        tmp_res = stringr::str_replace_all(tmp_res, 's', '/3600+')
        tmp_res = stringr::str_replace_all(tmp_res, 'd', '*24+')
        tmp_res = stringr::str_replace_all(tmp_res, 'm', '/60+')
        tmp_res = stringr::str_replace_all(tmp_res, 'h', '+')
        tmp_res = stringr::str_replace_all(tmp_res, ' ', '')
        tmp_res = stringr::str_replace_all(tmp_res, "\\+$", '')
        tmp_res = stringr::str_replace_all(tmp_res, "GB", "")
        tmp_res = eval(parse(text=tmp_res))
        return(tmp_res)
    }

# FIND ALL SLURM PROCESSES
    all_slurms = system('ls /project/holstegelab/Software/snakemake_pipeline/slurms_outputs/*out', intern = T)

# 1. DCACHE COPY (only runtime)
    # dcache copy
    dcache_copy = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep dcache"), intern = T)
        if (length(tmp_res) >0){
            tmp_res = system(paste0("tail ", f, " | grep 'Elapsed time'"), intern = T)
            if (length(tmp_res) > 0){
                tmp_res = stringr::str_split_fixed(tmp_res, ':', 2)[, 2]
                # make the time numeric
                tmp_res = transl(tmp_res)
                dcache_copy = rbind(dcache_copy, data.frame(id = slurm_id, runtime = tmp_res))
            }
        }
    }
    mean_runtime_dcache = mean(dcache_copy$runtime)
    median_runtime_dcache = median(dcache_copy$runtime)

# 2. CCS ALGORITHM
    all_ccs_log = system(paste0('find /project/holstegelab/Share/pacbio/data_processed -name "*.ccs.log"'), intern = T)
    # grep run time, cpu time and memory usage
    ccs_log = data.frame()
    for (f in all_ccs_log){
        cat(paste0('** processing ', f, '                                                                \r'))
        # runtime
        runtime = stringr::str_split_fixed(stringr::str_split_fixed(system(paste0('grep "Run Time" ', f), intern = T), '\\|', 4)[, 4], ':', 2)[, 2]
        if (length(runtime) >0){
            runtime = stringr::str_replace_all(runtime, ' ', '')
            # cpu time
            cputime = stringr::str_split_fixed(stringr::str_split_fixed(system(paste0('grep "CPU Time" ', f), intern = T), '\\|', 4)[, 4], ':', 2)[, 2]
            cputime = stringr::str_replace_all(cputime, ' ', '')
            # memory
            memory = stringr::str_split_fixed(stringr::str_split_fixed(system(paste0('grep "RSS" ', f), intern = T), '\\|', 4)[, 4], ':', 2)[, 2]
            memory = stringr::str_replace_all(memory, ' ', '')
            # save
            ccs_log = rbind(ccs_log, data.frame(file = f, runtime = runtime, cputime = cputime, memory = memory))
        }
    }
    # need to change the time -- convert d into days, h into hours, ms into minutes
    # runtime
    ccs_log$runtime_numeric = ccs_log$runtime
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 'ms', '/3600/1000+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 'us', '/3600/1000000+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 'ns', '/3600/1000000000+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 'm', '/60+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 'd', '*24+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 'h', '+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, 's', '/3600+')
    ccs_log$runtime_numeric = stringr::str_replace_all(ccs_log$runtime_numeric, "\\+$", "")
    ccs_log$runtime_numeric_hours = NULL
    for (i in 1:nrow(ccs_log)){
        ccs_log$runtime_numeric_hours[i] = eval(parse(text=ccs_log$runtime_numeric[i]))
    }
    # cpu time
    ccs_log$cputime_numeric = ccs_log$cputime
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 'ms', '/3600/1000+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 'us', '/3600/1000000+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 'ns', '/3600/1000000000+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 'm', '/60+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 'd', '*24+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 'h', '+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, 's', '/3600+')
    ccs_log$cputime_numeric = stringr::str_replace_all(ccs_log$cputime_numeric, "\\+$", "")
    ccs_log$cputime_numeric_hours = NULL
    for (i in 1:nrow(ccs_log)){
        ccs_log$cputime_numeric_hours[i] = eval(parse(text=ccs_log$cputime_numeric[i]))
    }
    # memory
    ccs_log$memory_gb = stringr::str_replace_all(ccs_log$memory, 'GB', '')
    # plot
    mean_runtime_ccs = mean(ccs_log$runtime_numeric_hours)
    mean_cputime_ccs = mean(ccs_log$cputime_numeric_hours)
    mean_memory_ccs = mean(as.numeric(ccs_log$memory_gb))
    # reduced representation
    ccs_log_red = ccs_log[, c('file', 'runtime_numeric_hours', 'cputime_numeric_hours', 'memory_gb')]
    colnames(ccs_log_red) = c('id', 'runtime', 'cputime', 'memory')

# 3. PRIMROSE (only runtime)
    primrose = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w primrose"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "idx"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            primrose = rbind(primrose, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(primrose$runtime))
    median(na.omit(primrose$runtime))

# 4. EXTRACT-CCS (only runtime)
    extract_ccs_nonccs = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w extract_ccs_nonccs"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "idx"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            extract_ccs_nonccs = rbind(extract_ccs_nonccs, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(extract_ccs_nonccs$runtime))
    median(na.omit(extract_ccs_nonccs$runtime))

# 5. CCS ALIGNMENT -- HG38
    # dcache copy
    ccs_aln = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w align_ccs"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':')))]]
            tmp_cputime = unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':')))]]
            tmp_memory = unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':')))]]
            if (length(tmp_runtime) > 0){
                # make the times numeric
                runtime_h = transl(tmp_runtime)
                cputime_h = transl(tmp_cputime)
                memory_gb = transl(tmp_memory)
                ccs_aln = rbind(ccs_aln, data.frame(id = slurm_id, runtime = runtime_h, cputime = cputime_h, memory = memory_gb))
            }
        }
    }
    mean_runtime_ccs_aln = mean(ccs_aln$runtime)
    mean_cputime_ccs_aln = mean(ccs_aln$cputime)
    mean_memory_ccs_aln = mean(ccs_aln$memory)
    median_runtime_ccs_aln = median(ccs_aln$runtime)
    median_cputime_ccs_aln = median(ccs_aln$cputime)
    median_memory_ccs_aln = median(ccs_aln$memory)
    
# 6. CCS ALIGNMENT -- CHM13
    # dcache copy
    ccs_aln_chm13 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep align_ccs_chm13"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':')))]]
            tmp_cputime = unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':')))]]
            tmp_memory = unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':')))]]
            if (length(tmp_runtime) > 0){
                # make the times numeric
                runtime_h = transl(tmp_runtime)
                cputime_h = transl(tmp_cputime)
                memory_gb = transl(tmp_memory)
                ccs_aln_chm13 = rbind(ccs_aln_chm13, data.frame(id = slurm_id, runtime = runtime_h, cputime = cputime_h, memory = memory_gb))
            }
        }
    }
    mean_runtime_ccs_aln_chm13 = mean(ccs_aln_chm13$runtime)
    mean_cputime_ccs_aln_chm13 = mean(ccs_aln_chm13$cputime)
    mean_memory_ccs_aln_chm13 = mean(ccs_aln_chm13$memory)
    median_runtime_ccs_aln_chm13 = median(ccs_aln_chm13$runtime)
    median_cputime_ccs_aln_chm13 = median(ccs_aln_chm13$cputime)
    median_memory_ccs_aln_chm13 = median(ccs_aln_chm13$memory)
    
# 7. NON-CCS ALIGNMENT -- HG38
    # dcache copy
    nonccs_aln = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w align_nonccs"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':')))]]
            tmp_cputime = unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':')))]]
            tmp_memory = unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':')))]]
            if (length(tmp_runtime) > 0){
                # make the times numeric
                runtime_h = transl(tmp_runtime)
                cputime_h = transl(tmp_cputime)
                memory_gb = transl(tmp_memory)
                nonccs_aln = rbind(nonccs_aln, data.frame(id = slurm_id, runtime = runtime_h, cputime = cputime_h, memory = memory_gb))
            }
        }
    }
    mean_runtime_nonccs_aln = mean(nonccs_aln$runtime)
    mean_cputime_nonccs_aln = mean(nonccs_aln$cputime)
    mean_memory_nonccs_aln = mean(nonccs_aln$memory)
    median_runtime_nonccs_aln = median(nonccs_aln$runtime)
    median_cputime_nonccs_aln = median(nonccs_aln$cputime)
    median_memory_nonccs_aln = median(nonccs_aln$memory)

# 8. NON-CCS ALIGNMENT -- CHM13
    # dcache copy
    nonccs_aln_chm13 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w align_nonccs_chm13"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Run Time'"), intern = T), ':')))]]
            tmp_cputime = unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'CPU Time'"), intern = T), ':')))]]
            tmp_memory = unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':'))[[length(unlist(strsplit(system(paste0("tail ", f, " | grep 'Peak RSS'"), intern = T), ':')))]]
            if (length(tmp_runtime) > 0){
                # make the times numeric
                runtime_h = transl(tmp_runtime)
                cputime_h = transl(tmp_cputime)
                memory_gb = transl(tmp_memory)
                nonccs_aln_chm13 = rbind(nonccs_aln_chm13, data.frame(id = slurm_id, runtime = runtime_h, cputime = cputime_h, memory = memory_gb))
            }
        }
    }
    mean_runtime_nonccs_aln_chm13 = mean(nonccs_aln_chm13$runtime)
    mean_cputime_nonccs_aln_chm13 = mean(nonccs_aln_chm13$cputime)
    mean_memory_nonccs_aln_chm13 = mean(nonccs_aln_chm13$memory)
    median_runtime_nonccs_aln_chm13 = median(nonccs_aln_chm13$runtime)
    median_cputime_nonccs_aln_chm13 = median(nonccs_aln_chm13$cputime)
    median_memory_nonccs_aln_chm13 = median(nonccs_aln_chm13$memory)

# 9. SAMPLE CHECK (only runtime)
    sample_check = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w sample_check"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "idx"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            sample_check = rbind(sample_check, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(sample_check$runtime))
    median(na.omit(sample_check$runtime))

# 10. COVERAGE STAT (only runtime)
    coverage_summary = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w coverage_summary"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "idx"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            coverage_summary = rbind(coverage_summary, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(coverage_summary$runtime))
    median(na.omit(coverage_summary$runtime))

# 11. MERGE CCS -- HG38 (only runtime)
    merge_hifi_hg38 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w merge_hifi_hg38"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "prefix"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            merge_hifi_hg38 = rbind(merge_hifi_hg38, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(merge_hifi_hg38$runtime))
    median(na.omit(merge_hifi_hg38$runtime))

# 12. MERGE CCS -- CHM13 (only runtime)
    merge_hifi_chm13 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w merge_hifi_chm13"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "prefix"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            merge_hifi_chm13 = rbind(merge_hifi_chm13, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(merge_hifi_chm13$runtime))
    median(na.omit(merge_hifi_chm13$runtime))

# 13. MERGE NON-CCS -- HG38 (only runtime)
    merge_nonhifi_hg38 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w merge_nonhifi_hg38"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "prefix"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            merge_nonhifi_hg38 = rbind(merge_nonhifi_hg38, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(merge_nonhifi_hg38$runtime))
    median(na.omit(merge_nonhifi_hg38$runtime))

# 14. MERGE NON-CCS -- CHM13 (only runtime)
    merge_nonhifi_chm13 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w merge_nonhifi_chm13"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "prefix"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            merge_nonhifi_chm13 = rbind(merge_nonhifi_chm13, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(merge_nonhifi_chm13$runtime))
    median(na.omit(merge_nonhifi_chm13$runtime))

# 15. PILEUP -- HG38 (only runtime)
    pileup_analysis_hg38 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w pileup_analysis_hg38"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "prefix" | grep -v "return" | grep -v "node" | grep -v "?"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            pileup_analysis_hg38 = rbind(pileup_analysis_hg38, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(pileup_analysis_hg38$runtime))
    median(na.omit(pileup_analysis_hg38$runtime))

# 16. PILEUP -- CHM13 (only runtime)
    pileup_analysis_chm13 = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w pileup_analysis_chm13"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "\\[" ', f, ' | grep -v "prefix" | grep -v "return" | grep -v "node" | grep -v "?"'), intern = T), ' ', 5)[, 4]
            diff_hours = as.numeric(stringr::str_split_fixed(difftime(strptime(tmp_runtime[2], "%H:%M:%S"), strptime(tmp_runtime[1], "%H:%M:%S"), units = 'mins'), ' ', 5)[, 1])/60
            if (!is.na(diff_hours) & diff_hours <0) { diff_hours = 24 - abs(diff_hours) }
            pileup_analysis_chm13 = rbind(pileup_analysis_chm13, data.frame(id = slurm_id, runtime = diff_hours))
        }
    }
    mean(na.omit(pileup_analysis_chm13$runtime))
    median(na.omit(pileup_analysis_chm13$runtime))

# 17. DEEPVARIANT 
    deepvariant = data.frame()
    for (f in all_slurms){
        slurm_id = unlist(strsplit(f, '/'))[[length(unlist(strsplit(f, '/')))]]
        cat(paste0('** processing ', slurm_id, '                   \r'))
        tmp_res = system(paste0('head ', f, " | grep -w deepvariant_hg38"), intern = T)
        if (length(tmp_res) >0){
            tmp_runtime = stringr::str_split_fixed(system(paste0('grep "real" ', f, ' | grep -v "time"'), intern = T), '\t', 2)[, 2]
            runtimes = c(); for (x in tmp_runtime){ runtimes = c(runtimes, transl(x)) }; runtime_sum = sum(runtimes)
            tmp_cputime = stringr::str_split_fixed(system(paste0('grep "user" ', f, ' | grep -v "time"'), intern = T), '\t', 2)[, 2]
            tmp_cputime_2 = stringr::str_split_fixed(system(paste0('grep "sys" ', f, ' | grep -v "time" | grep -v "support"'), intern = T), '\t', 2)[, 2]
            cputimes = c(); for (x in c(tmp_cputime, tmp_cputime_2)){ cputimes = c(cputimes, transl(x)) }; cputimes_sum = sum(cputimes)
            deepvariant = rbind(deepvariant, data.frame(id = slurm_id, runtime = runtime_sum, cputime = cputimes_sum))
        }
    }
    mean(na.omit(deepvariant$runtime))
    median(na.omit(deepvariant$runtime))
    mean(na.omit(deepvariant$cputime))
    median(na.omit(deepvariant$cputime))
###############################################
# 18. COMBINE DATA
    # remove NAs before
    dcache_copy = dcache_copy[!is.na(dcache_copy$runtime),]
    ccs_log_red = ccs_log_red[!is.na(ccs_log_red$runtime),]
    primrose = primrose[!is.na(primrose$runtime),]
    extract_ccs_nonccs = extract_ccs_nonccs[!is.na(extract_ccs_nonccs$runtime),]
    ccs_aln = ccs_aln[!is.na(ccs_aln$runtime),]
    ccs_aln_chm13 = ccs_aln_chm13[!is.na(ccs_aln_chm13$runtime),]
    sample_check = sample_check[!is.na(sample_check$runtime),]
    coverage_summary = coverage_summary[!is.na(coverage_summary$runtime),]
    merge_hifi_hg38 = merge_hifi_hg38[!is.na(merge_hifi_hg38$runtime),]
    merge_hifi_chm13 = merge_hifi_chm13[!is.na(merge_hifi_chm13$runtime),]
    merge_nonhifi_hg38 = merge_nonhifi_hg38[!is.na(merge_nonhifi_hg38$runtime),]
    merge_nonhifi_chm13 = merge_nonhifi_chm13[!is.na(merge_nonhifi_chm13$runtime),]
    pileup_analysis_hg38 = pileup_analysis_hg38[!is.na(pileup_analysis_hg38$runtime),]
    pileup_analysis_chm13 = pileup_analysis_chm13[!is.na(pileup_analysis_chm13$runtime),]
    deepvariant = deepvariant[!is.na(deepvariant$runtime),]
    # make sure the column names are the same for the rbind
    dcache_copy$operation = 'dcache copy'
    ccs_log_red$operation = 'ccs'
    primrose$operation = 'primrose'
    extract_ccs_nonccs$operation = 'extract ccs'
    ccs_aln$operation = 'ccs alignment hg38'
    ccs_aln_chm13$operation = 'ccs alignment chm13'
    nonccs_aln$operation = 'nonccs alignment hg38'
    nonccs_aln_chm13$operation = 'nonccs alignment chm13'
    sample_check$operation = 'sample check'
    coverage_summary$operation = 'coverage summary'
    merge_hifi_hg38$operation = 'merge ccs hg38'
    merge_hifi_chm13$operation = 'merge ccs chm13'
    merge_nonhifi_hg38$operation = 'merge nonccs hg38'
    merge_nonhifi_chm13$operation = 'merge nonccs chm13'
    pileup_analysis_hg38$operation = 'pileup hg38'
    pileup_analysis_chm13$operation = 'pileup chm13'
    deepvariant$operation = 'deepvariant'
    combined_data = rbind.fill(dcache_copy, ccs_log_red, primrose, extract_ccs_nonccs, ccs_aln, ccs_aln_chm13, nonccs_aln, nonccs_aln_chm13, sample_check, coverage_summary, merge_hifi_hg38, merge_hifi_chm13, merge_nonhifi_hg38, merge_nonhifi_chm13, pileup_analysis_hg38, pileup_analysis_chm13, deepvariant)
    combined_data$operation = factor(combined_data$operation, levels = c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 
        'nonccs alignment chm13', 'sample check', 'coverage summary', 'merge ccs hg38', 'merge ccs chm13', 'merge nonccs hg38', 'merge nonccs chm13', 'pileup hg38', 'pileup chm13', 'deepvariant'))
    combined_data$memory = as.numeric(combined_data$memory)
    # plot -- runtime
    runtime_p = ggplot(combined_data, aes(x = operation, y = log10(runtime), fill = operation)) + geom_boxplot() + xlab('Operations') + ylab('log10 (Runtime in Hours)') + ggtitle('Runtime across operations') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    # plot -- cputime
    cputime_p = ggplot(combined_data, aes(x = operation, y = log10(cputime), fill = operation)) + geom_boxplot() + xlab('Operations') + ylab('log10 (CPU Hours)') + ggtitle('CPU time across operations') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    # plot -- memory
    memory_p = ggplot(combined_data, aes(x = operation, y = memory, fill = operation)) + geom_boxplot() + xlab('Operations') + ylab('Memory (GB)') + ggtitle('Memory across operations') + 
        theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    # combine plots
    pdf('plot_processing_runtime_cputime_memory.pdf', height = 12, width = 8)
    combined_plot = egg::ggarrange(runtime_p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
        cputime_p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()),
        memory_p, nrow = 3, ncol = 1)
    dev.off()
# 19. ESTIMATE RUN TIME AND RESOURCES PER RUN
    # runtime
    groups_mean_runtime = aggregate(combined_data$runtime, list(combined_data$operation), mean); colnames(groups_mean_runtime) = c('variable', 'mean')
    groups_median_runtime = aggregate(combined_data$runtime, list(combined_data$operation), median); colnames(groups_median_runtime) = c('variable', 'median')
    # cputime
    groups_mean_cputime = aggregate(combined_data$cputime, list(combined_data$operation), mean); colnames(groups_mean_cputime) = c('variable', 'mean')
    groups_median_cputime = aggregate(combined_data$cputime, list(combined_data$operation), median); colnames(groups_median_cputime) = c('variable', 'median')
    # memory
    groups_mean_memory = aggregate(combined_data$memory, list(combined_data$operation), mean); colnames(groups_mean_memory) = c('variable', 'mean')
    groups_median_memory = aggregate(combined_data$memory, list(combined_data$operation), median); colnames(groups_median_memory) = c('variable', 'median')
    # a single run is based on all steps until the merging
    # runtime
    single_sample_runtime_mean = sum(groups_mean_runtime$mean[which(groups_mean_runtime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary'))])
    single_sample_runtime_median = sum(groups_median_runtime$median[which(groups_median_runtime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary'))])
    # cputime
    single_sample_cputime_mean = sum(na.omit(groups_mean_cputime$mean[which(groups_mean_cputime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary'))]))
    single_sample_cputime_median = sum(na.omit(groups_median_cputime$median[which(groups_median_cputime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary'))]))
# 20. For a full sample we can assume 2 smrt cells + merging and so on
    # runtime
    smrt_cell_sample_mean_runtime = 2*single_sample_runtime_mean
    smrt_cell_sample_median_runtime = 2*single_sample_runtime_median
    # cputime
    smrt_cell_sample_mean_cputime = 2*single_sample_cputime_mean
    smrt_cell_sample_median_cputime = 2*single_sample_cputime_median
    # add the merging and deepvariant
    # runtime
    single_sample_runtime_mean_merge = sum(groups_mean_runtime$mean[which(!(groups_mean_runtime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary')))])
    single_sample_runtime_median_merge = sum(groups_median_runtime$median[which(!(groups_median_runtime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary')))])
    # cputime
    single_sample_cputime_mean_merge = sum(na.omit(groups_mean_cputime$mean[which(!(groups_mean_cputime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary')))]))
    single_sample_cputime_median_merge = sum(na.omit(groups_median_cputime$median[which(!(groups_median_cputime$variable %in% c('dcache copy', 'ccs', 'primrose', 'extract ccs', 'ccs alignment hg38', 'ccs alignment chm13', 'nonccs alignment hg38', 'nonccs alignment chm13', 'sample check', 'coverage summary')))]))
    # finally we do the sum
    # runtime
    total_time_per_sample_mean = smrt_cell_sample_mean_runtime + single_sample_runtime_mean_merge
    total_time_per_sample_median = smrt_cell_sample_median_runtime + single_sample_runtime_median_merge
    # cputime
    total_cpu_per_sample_mean = smrt_cell_sample_mean_cputime + single_sample_cputime_mean_merge
    total_cpu_per_sample_median = smrt_cell_sample_median_cputime + single_sample_cputime_median_merge
    # plot this
    df = data.frame(time = c(single_sample_runtime_mean, single_sample_cputime_mean, total_time_per_sample_mean, total_cpu_per_sample_mean), 
                    type = c('Hours', 'CPU Hours', 'Hours', 'CPU Hours'), per = c('SMRT cell run', 'SMRT cell run', 'Sample (2x SMRT runs)', 'Sample (2x SMRT runs)'))
    p = ggplot(df, aes(x = type, y = time, fill = type)) + geom_bar(stat = 'identity') + xlab('') + ylab('Time (Hours)') + ggtitle('Estimated Hours at run and sample level')
    p = p + facet_grid(cols = vars(per)) + theme(legend.position = 'none') + geom_hline(yintercept=100, linetype="dashed", color = "red") + geom_hline(yintercept=50, linetype="dashed", color = "blue")
    pdf('plot_summary_per_sample.pdf', height = 8, width = 12)
        print(p)
    dev.off()

# 21. Evaluate storage
    # list all processed files
    all_proces_bbc = system('find /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/ -name "m*"', intern = T)
    all_proces_adchc = system('find /project/holstegelab/Share/pacbio/data_processed/blood_brain_child/ -name "m*"', intern = T)
    all_proces_nijmegen = system('find /project/holstegelab/Share/pacbio/data_processed/nijmegen/ -name "m*"', intern = T)
    all_proces = c(all_proces_bbc, all_proces_adchc, all_proces_nijmegen)
    # clean this data by keeping only the filename
    all_smrt = c()
    tmp = strsplit(all_proces, '/')
    for (x in tmp){ 
        smrt_id = stringr::str_split_fixed(x[length(x)], '\\.', 2)[, 1]
        all_smrt = c(all_smrt, smrt_id)
    }
    all_smrt = all_smrt[!duplicated(all_smrt)]
    # estimate size per smrt run
    storage = data.frame()
    for (x in all_smrt){
        cat(paste0('** working on ', x, '                    \r'))
        flist = system(paste0('find /project/holstegelab/Share/pacbio/data_processed/ -name "', x, '.ccs*"'), intern = T)
        sizes = c()
        for (f in flist){
            size = system(paste0("ls -lh ", f, " | sed 's/ /\t/g' | cut -f5,5"), intern = T)
            if (length(grep('G', size)) >0){ 
                size = as.numeric(stringr::str_replace_all(size, 'G', '')) 
            } else if (length(grep('M', size)) >0){ 
                size = as.numeric(stringr::str_replace_all(size, 'M', ''))/1000 
            } else if (length(grep('K', size)) >0){ 
                size = as.numeric(stringr::str_replace_all(size, 'K', ''))/1000000 
            } else {
                size = as.numeric(size)/1000000000 
            }
            sizes = c(sizes, size)
        }
        files_size = sum(sizes)
        storage = rbind(storage, data.frame(id = x, total = files_size))
    }
    mean(storage$total)
    median(storage$total)
# 22. save workspace
    save.image('20221007_workspace_estimation.RData')









