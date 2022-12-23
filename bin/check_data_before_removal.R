# SCRIPT TO CHECK WHETHER DATA CAN BE REMOVED FROM CLINICAL GENETICS SERVER
# TO BE REMOVED, DATA NEEDS TO BE:
# 1. ON DCACHE
# 2. DATA SHOULD ALSO BE PROCESSED

# path to data from KG to be removed (this should be a list of runs)
kg_path = 'smrt_cells_to_be_remove_KG_20221223.txt'
kg = read.table(kg_path, h=F, stringsAsFactors=F)

# for each file, check if a link exists in dcache
run_stat = data.frame()
for (smrt in kg$V1){
    if (file.exists(paste0('/home/holstegelab-ntesi/dcache/tape/pacbio/', smrt)) == TRUE){
        # if file exists, check if file was correctly processed
        smrt_id = system(paste0("ls /home/holstegelab-ntesi/dcache/tape/pacbio/", smrt, "/*transferdone"), intern = T)
        smrt_id = stringr::str_replace_all(unlist(strsplit(smrt_id, '/'))[[length(unlist(strsplit(smrt_id, '/')))]], '\\.transferdone', '')
        processed_files = system(paste0('ls /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/', smrt_id, '.*'), intern = T)
        processed_files_dcache = system(paste0('ls ~/dcache/tape/data_processed/ccs/ad_centenarians/', smrt_id, '.*'), intern = T)
        processed_files_all = c(processed_files, processed_files_dcache)
        # save stats
        if (length(processed_files_all) >0){
            tmp_stat = data.frame(kg_id = smrt, smrt_id = smrt_id, on_dcache = TRUE, processed = TRUE, n_downstream_files = length(processed_files))
        } else {
            tmp_stat = data.frame(kg_id = smrt, smrt_id = smrt_id, on_dcache = TRUE, processed = FALSE, n_downstream_files = 0)
        }
    } else {
        tmp_stat = data.frame(kg_id = smrt, smrt_id = smrt_id, on_dcache = FALSE, processed = FALSE, n_downstream_files = 0)
    }
    run_stat = rbind(run_stat, tmp_stat)
}
# ok this worked and for most data is ok
table(run_stat$on_dcache, exclude=F) # --> if this is = to the number of runs to check, we are ok, all runs are on dcache
table(run_stat$processed, run_stat$n_downstream_files)
# some runs have not been processed -- maybe yield was too low?
yield_missing = data.frame()
for (i in 1:nrow(run_stat)){
    if (run_stat$processed[i] == FALSE){
        yield = system(paste0("ls -lh /home/holstegelab-ntesi/dcache/tape/pacbio/", run_stat$kg_id[i], "/*subreads.bam | sed 's/ /\t/g' | cut -f5,5"), intern = T)
        if (length(yield) == 0) yield = NA
        yield_missing = rbind(yield_missing, data.frame(smrt_kg = run_stat$kg_id[i], yield = yield))
    }
}
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


