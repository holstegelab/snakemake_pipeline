# SCRIPT TO CHECK WHETHER DATA CAN BE REMOVED FROM CLINICAL GENETICS SERVER
# TO BE REMOVED, DATA NEEDS TO BE:
# 1. ON DCACHE
# 2. DATA SHOULD ALSO BE PROCESSED

# path to data from KG to be removed (this should be a list of runs)
kg_path = 'smrt_cells_to_be_remove_KG_20230720.txt'
kg = read.table(kg_path, h=F, stringsAsFactors=F)

# for each file, check if a link exists in dcache
run_stat = data.frame()
for (smrt in kg$V1){
    print(paste0('** working on run --> ', smrt))
    if (file.exists(paste0('/home/holstegelab-ntesi/dcache/tape/pacbio/', smrt)) == TRUE){
        # if file exists, check if file was correctly processed
        smrt_id = system(paste0("ls /home/holstegelab-ntesi/dcache/tape/pacbio/", smrt, "/*transferdone"), intern = T)
        if (length(smrt_id) > 0){
            smrt_id = stringr::str_replace_all(unlist(strsplit(smrt_id, '/'))[[length(unlist(strsplit(smrt_id, '/')))]], '\\.transferdone', '')
            processed_files1 = system(paste0('ls /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/', smrt_id, '.*'), intern = T)
            processed_files2 = system(paste0('ls /project/holstegelab/Share/pacbio/data_processed/Anke_samples/', smrt_id, '.*'), intern = T)
            processed_files3 = system(paste0('ls /project/holstegelab/Share/pacbio/data_processed/other_samples/', smrt_id, '.*'), intern = T)
            processed_files_dcache1 = system(paste0('ls ~/dcache_processed/ccs/ad_centenarians/', smrt_id, '.*'), intern = T)
            processed_files_dcache2 = system(paste0('ls ~/dcache_processed/ccs/other/', smrt_id, '.*'), intern = T)
            processed_files_dcache3 = system(paste0('ls ~/dcache_processed/ccs/blood_brain_child/', smrt_id, '.*'), intern = T)
            processed_files_all = c(processed_files1, processed_files2, processed_files3, processed_files_dcache1, processed_files_dcache2, processed_files_dcache3)
            # save stats
            if (length(processed_files_all) >0){
                tmp_stat = data.frame(kg_id = smrt, smrt_id = smrt_id, on_dcache = TRUE, processed = TRUE, n_downstream_files = length(processed_files_all))
            } else {
                tmp_stat = data.frame(kg_id = smrt, smrt_id = smrt_id, on_dcache = TRUE, processed = FALSE, n_downstream_files = 0)
            }
        } else {
            tmp_stat = data.frame(kg_id = smrt, smrt_id = NA, on_dcache = FALSE, processed = FALSE, n_downstream_files = 0)    
        }
    } else {
        tmp_stat = data.frame(kg_id = smrt, smrt_id = smrt_id, on_dcache = FALSE, processed = FALSE, n_downstream_files = 0)
    }
    run_stat = rbind(run_stat, tmp_stat)
}
# ok this worked and for most data is ok
table(run_stat$on_dcache) # --> if this is = to the number of runs to check, we are ok, all runs are on dcache
table(run_stat$processed, run_stat$n_downstream_files)
# show runs not on dcache
run_stat[which(run_stat$on_dcache == FALSE),]
run_stat[which(run_stat$processed == FALSE),]
# some runs have not been processed -- maybe yield was too low?
yield_missing = data.frame()
for (i in 1:nrow(run_stat)){
    yield = system(paste0("ls -lh ~/dcache/tape/pacbio/", run_stat$kg_id[i], "/*subreads.bam | sed 's/ /\t/g' | cut -f5,5"), intern = T)
    if (length(yield) == 0) yield = NA
    yield_missing = rbind(yield_missing, data.frame(smrt_kg = run_stat$kg_id[i], yield = yield))
}
summary_runs = merge(run_stat, yield_missing, by.x = 'kg_id', by.y = 'smrt_kg')
summary_runs[which(summary_runs$processed == FALSE),]
# r64050e_20220729_121439/1_A01 --> 560G --> Control DNA
# r64367e_20220520_141902/2_B01 --> NA --> no raw subreads -- only ccs
# r64367e_20220513_112845/3_C01 --> 286G --> processed but it's Anke's sample
# r64367e_20220812_115756/1_A01 --> 953G --> Control DNA
# r64050e_20220812_115528/1_A01 --> 1.1T --> Control DNA
# r64367e_20220729_121556/1_A01 --> 906G --> Control DNA
#########################################################
# update of 2022-12-23
# all good: low quality smrt cells have not been processed, rest is ok
# everything is correctly on dcache
#########################################################
#########################################################
# update of 2023-02-08
# all good: low quality smrt cells have not been processed, rest is ok
# everything is correctly on dcache
#########################################################
# update of 2023-03-20
# of 13 SMRT cells, all of them are on dcache
# 8 were processed and data is ready
# 5 were not processed as the subreads.bam was basically empty. Check with Bilgehan whether this was an issue like in the past.
#########################################################
# update of 2023-06-26
# of 65 SMRT cells, all of them are on dcache
# all the good ones have been processed
#########################################################
# update of 2023-07-20
# of 40 SMRT cells, all of them are on dcache
# all the good ones have been processed

