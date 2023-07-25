####################################
####################################
# MAKE AN OVERVIEW OF THE SEQUENCING STATUS
    # LIBRARIES
    library(data.table)
    library(plyr)
    library(openxlsx)
    library(stringr)
    library(ggplot2)

    # READ ALL SAMPLES WE WERE SUPPOSED TO DO
    all_chc = fread('chc_for_pacbio_withSCORE_and_KGnumb.txt', h=T)
    all_ad = fread('20200724_final_list_ADcases_completed.txt', h=T)
    # COMBINE THEM
    colnames(all_chc) = c('ID_100plus', 'ID_GWAS', 'pheno', 'gender', 'age', 'pca_outliers_gwas', 'pbmc', 'children', 'children_partners', 'siblings', 'trios', 'brain', 'score', 'KG-number')
    colnames(all_ad) = c('I_ID', 'pheno', 'info', 'gender', 'brain_year', 'apoe', 'age_death', 'age_diag', 'last_age_available', 'brain', 'KG-number')
    # for AD cases, need to add the GWAS ID
    load('/project/holstegelab/Share/gwas_array/mapping_files/phenotypes_20211027_All.Rdata')
    tmp_pheno = pheno_final_raw[, c('I_ID', 'ID_GWAS')]
    all_ad = merge(all_ad, tmp_pheno, by = 'I_ID', all.x = T)
    # there are 3 duplicates: that is, 3 samples (AD) with 2 ID_GWAS. 2 of them make sense, as they were genotyped multiple times, 1 is more controversial as it's from the same genotyping run but different phenotypes (male/female, I_ID = 2100). Keep all duplicates in for now.
    # merge    
    all_samples_expected = rbind.fill(all_chc, all_ad)
    table(all_samples_expected$pheno)
    # we started of with 353 AD cases (3 duplicates) and 356 centenarians

    # READ ALL PLATES
    plates = data.frame()
    # plate 1 is in a different data, do it separately
        plate1 = read.xlsx('Plate1_QC.xlsx', sheet = 1)
        # restrict to plate, biobank, barcode and KG-number
        plate1 = plate1[, c('PlateIDtarget', 'SampleID', 'Barcode', 'Individual.ID')]
        # rename columns
        colnames(plate1) = c('plate', 'biobank', 'barcode', 'kg_number')
        # add plate id
        plate1$plate_id = 1
        # add to plates
        plates = rbind(plates, plate1)
    for (p in 2:7){
        # read sheet
        tmp = read.xlsx('All_plates_QC.xlsx', sheet = p)
        # restrict to plate, biobank, barcode and KG-number
        tmp = tmp[, c('Plate', 'Biobank.fractie', 'BC', 'KG-number')]
        # rename columns
        colnames(tmp) = c('plate', 'biobank', 'barcode', 'kg_number')
        # add plate id
        tmp$plate_id = p
        # add to plates
        plates = rbind(plates, tmp)
    }
    table(plates$plate_id)
    # the plates contain 653 samples in total

    # NOW READ ALL SEQUENCED SAMPLES BASED ON THE COVERAGE SUMMARY FILE
    merged_sequenced = fread('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', h=T, stringsAsFactors=F, sep="\t", dec = ',')
    merged_sequenced = merged_sequenced[!is.na(merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[!duplicated(merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[order(-merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[!duplicated(merged_sequenced$SAMPLE),]
    # also read blood-brain-child samples
    bbc_samples = fread('/project/holstegelab/Software/snakemake_pipeline/sample_check_data/20210615_overview_genetics.txt', h=T, stringsAsFactors=F)
    bbc_samples = bbc_samples[!is.na(bbc_samples$ID_PACBIO),]
    # modify here
    bbc_samples = bbc_samples[, c('ID_GWAS')]
    bbc_samples$DIAGNOSIS = 'Blood_Brain_Child'
    colnames(bbc_samples) = c('SAMPLE', 'DIAGNOSIS')
    merged_sequenced = rbind.fill(merged_sequenced, bbc_samples)
    table(merged_sequenced$DIAGNOSIS)
    # in total, ready for analysis there are 60 AD cases + 82 Centenarians + 10 more centenarians (blood-brain-child project) + 10 children (blood-brain-child project)

    # NOW LABEL THE SAMPLES SEQUENCED IN ALL_SAMPLES WE WERE SUPPOSED TO SEQUENCE
    table(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)
    # 152/162 are OK -- let's check the remaining 10
    non_matching = merged_sequenced[which(!(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)),]
    # add phenotypes
    non_matching_pheno = merge(non_matching, pheno_final_raw, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    # it's ok -- these 10 are the children of the centenarians
    # So all samples we sequenced until now were in the list of samples to be sequenced

    # NOW CHECK WITH RESPECT TO THE PLATES
    # add kg-number to the sequenced samples
    kg_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
    merged_sequenced = merge(merged_sequenced, kg_info, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    table(merged_sequenced$DIAGNOSIS)
    # merge sequenced samples with plate information
    merged_sequenced_plate = merge(merged_sequenced, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    # add phenotypes
    tmp_pheno = pheno_final_raw[, c('ID_GWAS', 'diagnosis', 'ID_100plus', 'I_ID')]
    merged_sequenced_plate = merge(merged_sequenced_plate, tmp_pheno, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    table(merged_sequenced_plate$DIAGNOSIS)
    # there are duplicates, we should clean them up
    dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
    merged_sequenced_plate = merged_sequenced_plate[order(-merged_sequenced_plate$COMBINED_COVERAGE),]
    merged_sequenced_plate = merged_sequenced_plate[!duplicated(merged_sequenced_plate$SAMPLE),] 
    dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
    # find samples wiith missing plate id
    merged_sequenced_plate[is.na(merged_sequenced_plate$plate_id),]
    table(merged_sequenced_plate$plate_id, exclude=F)
    merged_sequenced_plate$SMRT_CELL_N = factor(merged_sequenced_plate$SMRT_CELL_N)
    # 14 in total: 10 are centenarians from the blood-brain-child project, 4 are additional centenarians (not sure why they are not in the list)

    # PLOT THE COMBINED_COVERAGE OF SAMPLES
    ggplot(merged_sequenced_plate[which(merged_sequenced_plate$DIAGNOSIS != 'Blood_Brain_Child'),], aes(x = DIAGNOSIS, y = COMBINED_COVERAGE, fill = DIAGNOSIS)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = SMRT_CELL_N), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1")

    # LOOK AT SMRT CELLS >600GB YIELD
    smrt_coverage_p1 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-02-17_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    smrt_coverage_p2 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-03-10_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    # subset columns
    smrt_coverage_p2 = smrt_coverage_p2[, c('movie_id', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    smrt_coverage_p1 = smrt_coverage_p1[, c('SMRT_ID', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    colnames(smrt_coverage_p2)[1] = 'SMRT_ID'
    # combine
    smrt_coverage_all = rbind(smrt_coverage_p1, smrt_coverage_p2)
    # remove duplicates based on smrt id
    dups = smrt_coverage_all[duplicated(smrt_coverage_all$SMRT_ID),]
    smrt_coverage_all = smrt_coverage_all[!duplicated(smrt_coverage_all$SMRT_ID),]
    dim(smrt_coverage_all)
    # ok now exclude samples combined (from merged_sequenced_plate)
    smrt_coverage_all_notcombined = smrt_coverage_all[which(!(smrt_coverage_all$ID_GWAS %in% merged_sequenced_plate$SAMPLE)),]
    # also exclude samples from blood-brain-child project
    smrt_coverage_all_notcombined = smrt_coverage_all_notcombined[which(!(smrt_coverage_all_notcombined$ID_GWAS %in% bbc_samples$SAMPLE)),]
    table_occurrences = data.frame(table(smrt_coverage_all_notcombined$ID_GWAS))
    table_occurrences = table_occurrences[order(-table_occurrences$Freq),]
    # some things are off: for example, sample 1128423 (centenarian) was merged but there's no track of it
    # calculate total coverage of these
    table_occurrences$Cov = NA
    table_occurrences$Pheno = NA
    table_occurrences$Smrts = NA
    for (i in 1:nrow(table_occurrences)){
        tmp = smrt_coverage_all_notcombined[which(smrt_coverage_all_notcombined$ID_GWAS == table_occurrences$Var1[i]),]
        table_occurrences$Cov[i] = sum(tmp$GLOBAL_COVERAGE)
        table_occurrences$Pheno[i] = unique(tmp$diagnosis)
        table_occurrences$Smrts[i] = paste(tmp$SMRT_ID, collapse = ',')
    }
    table_occurrences = table_occurrences[order(-table_occurrences$Cov),]
    head(table_occurrences, 10)
    # get most updated list of samples merged
    samples_merged_chc = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/centenarian/*merged.hifi.chm13.bam', intern=T)
    samples_merged_ad = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/ad/*merged.hifi.chm13.bam', intern=T)
    df_chc = data.frame(sample = samples_merged_chc, pheno = 'chc')
    df_ad = data.frame(sample = samples_merged_ad, pheno = 'ad')
    df_combined = rbind(df_chc, df_ad)
    for (i in 1:nrow(df_combined)){
        tmp = str_replace_all(strsplit(df_combined$sample[i], '/')[[1]][length(strsplit(df_combined$sample[i], '/')[[1]])], '.merged.hifi.chm13.bam', '')
        df_combined$sample[i] = tmp
    }
    # exclude from the table of non combined any of these actually combined
    dim(table_occurrences)
    table_occurrences_filter = table_occurrences[which(!(table_occurrences$Var1 %in% df_combined$sample)), ]
    dim(table_occurrences_filter)
    # only 3 were done, good, probably need to update the merged files!
    head(table_occurrences_filter)
    # there are some things to keep in mind here:
    table_occurrences_filter$Reason = NA
    # 1. 11 runs (1, 2, 3, 4, 5, 7, 8) have to be merged because of the recently processed bam files from nijmegen
        table_occurrences_filter$Reason[c(1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 14)] = 'Recently processed Nijmegen'
    # 2. 2 run (6) have to be merged. this was done in amsterdam but was not uploaded until recently
        table_occurrences_filter$Reason[c(6, 10)] = 'Recently processed Amsterdam'
    # 3. 1 run can be removed (9) as this is a duplicated sample. The sample was processed but with the other ID.
        table_occurrences_filter = table_occurrences_filter[which(table_occurrences_filter$Var1 != '3196128'),]
    # then take out the samples with coverage >12
    coverage_enouth = table_occurrences_filter[which(table_occurrences_filter$Cov >= 12), ]
    dim(coverage_enouth)
    coverage_enouth$Freq = factor(coverage_enouth$Freq)
    # plot these
    ggplot(coverage_enouth[which(coverage_enouth$Pheno != 'Control_100plus'),], aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1")
    # also overlap with plates
    plates_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
    coverage_enouth = merge(coverage_enouth, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
    coverage_enouth_plate = merge(coverage_enouth, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    # check plates
    table(merged_sequenced_plate$plate_id, exclude=F)
    table(coverage_enouth_plate$plate_id, exclude=F)
    # also nice to plot the samples with not enough coverage
    not_enough = table_occurrences_filter[which(!(table_occurrences_filter$Var1 %in% coverage_enouth$Var1)),]
    not_enough$Pheno_short = NA; not_enough$Pheno_short[which(not_enough$Pheno %in% c('AD_othersite', 'Probable_AD'))] = 'Probable_AD'
    not_enough$Pheno_short[which(not_enough$Pheno %in% c('Centenarian', 'family_100plus'))] = 'Centenarian'
    not_enough$Pheno_short[is.na(not_enough$Pheno_short)] = 'Other'
    not_enough$Freq = factor(not_enough$Freq)

    ggplot(not_enough, aes(x = Pheno_short, y = Cov, fill = Pheno_short)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylim(0, 12)

    # next go to the Nijmegen stats
    nijm_logs = fread('Nijmegen_Log.txt', h=T, stringsAsFactors=F, dec=',')
    # exclude runa we have processed
    nijm_logs_missing = nijm_logs[which(!(nijm_logs$MovieID %in% smrt_coverage_all$SMRT_ID)),]
    # n=33 (makes sense, this is basically october and november)
    # plot yield and date
    nijm_logs_missing$date = str_split_fixed(nijm_logs_missing$MovieID, '_', 3)[, 2]
    nijm_logs_missing$date = paste0('20', nijm_logs_missing$date)
    nijm_logs_missing$date = as.Date(nijm_logs_missing$date, "%Y%m%d")
    colnames(nijm_logs_missing) = c('Smrt', 'Yield', 'Date')
    ggplot(nijm_logs_missing, aes(x=Date, y=Yield, color=Yield)) + geom_point(stat = 'identity', size=3) + ggtitle('Missing runs from Nijmegen only based on Log files') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
    # see how many runs with yield >600gb we have
    dim(nijm_logs_missing[which(nijm_logs_missing$Yield >= 600),])
    # maybe plot Yield over time (all smrt cells ever?)
    all_run_logs = fread('All_runs_Logs.txt', h=T, stringsAsFactors=F, dec=',')
    all_run_logs$date = str_split_fixed(all_run_logs$MovieID, '_', 3)[, 2]
    all_run_logs$date = paste0('20', all_run_logs$date)
    all_run_logs$date = as.Date(all_run_logs$date, "%Y%m%d")
    colnames(all_run_logs) = c('Smrt', 'Yield', 'Site', 'Date')
    all_run_logs = all_run_logs[!is.na(all_run_logs$Yield),]
    all_run_logs = all_run_logs[order(all_run_logs$Date),]
    ggplot(all_run_logs, aes(x=Date, y=Yield, color=Site)) + geom_point(stat = 'identity', size=3, alpha = 0.5) + geom_smooth(method = 'loess') + ggtitle('SMRT Yield over time and sequencing centers') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))

####################################
# 2023-03-16
# AIM IS TO GENERATE AN OVERVIEW AT THE SAMPLE LEVEL
    load('Desktop/20230313_samples_overview.RData')
    library(openxlsx)
    # This is based on the merged samples
    merged_sequenced
    # This is based on the single smrt cells
    table_occurrences_filter
    # check whether we expected to sequence these
    tofix = table_occurrences_filter[which(!(table_occurrences_filter$Var1 %in% all_samples_expected$ID_GWAS)),]
    table_occurrences_filter = table_occurrences_filter[which(!(table_occurrences_filter$Var1 %in% tofix$Var1)),]
    # something to be fixed
    updated_diagnosis = c('Centenarian', 'Probable_AD', 'NIID', 'Dementie_anders', 'family_100plus', 'NIID')
    updated_ids = c('8988312', '9665575', '20R1483', '4385663', '7797040', '20R1484')
    tofix$Var1 = updated_ids; tofix$Pheno = updated_diagnosis
    table_occurrences_filter = rbind(table_occurrences_filter, tofix)
    table_occurrences_filter_info = merge(table_occurrences_filter, pheno_final_raw, by.x = 'Var1', by.y = 'ID_GWAS')
    table_occurrences_filter_info2 = merge(table_occurrences_filter_info, all_samples_expected, by.x = 'Var1', by.y = 'ID_GWAS', all.x=T)
    
    # then we need to go to Nijmegen data
    new_nijm = read.xlsx('Desktop/Nijmegen_Sample_ID_and_RunIDs_Overview.xlsx')
    # exclude samples from the blood-brain-child project
    new_nijm = new_nijm[which(new_nijm$Project != 'BBC'),]
    # from this new list, exclude smrt cells already processed (in table occurrences)
    new_nijm$Done = NA
    for (i in 1:nrow(new_nijm)){
        # get run id
        tmp_smrt = new_nijm$MovieID[i]
        # check if was merged
        was_merged = merged_sequenced[grep(tmp_smrt, merged_sequenced$SMRT_CELLS),]
        was_proces = table_occurrences_filter[grep(tmp_smrt, table_occurrences_filter$Smrts),]
        if (nrow(was_merged) >0){
            new_nijm$Done[i] = 'merged'
        } else if (nrow(was_proces) >0){
            new_nijm$Done[i] = 'proces'
        } else{
            new_nijm$Done[i] = 'todo'
        }
    }
    table(new_nijm$Done, exclude=F)
    # then take only todo
    todo = new_nijm[which(new_nijm$Done == 'todo'),]
    todo$Sample = stringr::str_replace_all(todo$Sample, '-02', '')
    todo$Sample = stringr::str_replace_all(todo$Sample, '_19kb_30-07-2021', '')
    
    # read run stats
    run_stats = read.xlsx('Desktop/Copy_nijmegen_runs_from_okt.xlsx', sheet = 2)
    table(todo$MovieID %in% run_stats$context)
    not_matching = todo[which(!todo$MovieID %in% run_stats$context),]
    todo = merge(todo, run_stats, by.x = 'MovieID', by.y = 'context', all.x=T)
    # 16 do not match between the run statistics from Nijmegen: this can be because the run statistics (run_stats) are only from october onwards, while the other file (todo) also includes other months, so that's ok
    table(todo$Sample %in% plates$barcode)
    # 7 samples sequenced were not in the plates -- maybe some change we did afterwards?
    mismatch = todo[which(!(todo$Sample %in% plates$barcode)),]
        # solve these 7
        missing1 = c('KG-011524', '17R1457',	'FR34137283', '202204011P1V', 'D01')
        missing2 = c('ND-001200', '11D01145', 'FR34137218', '202204011P1V', 'A01')
        missing3 = c('ND-000882', '10D04356', 'FR34137285', '202204011P1V', 'G04')
        missing4 = c('ND-000025', '09D04551', 'FR34137227', '202204011P1V', 'B06')
        missing5 = c('KG-011559', '17R1369', 'FR34137302', '202204011P1V', 'H01')
        df_missing = data.frame(plate = c(missing1[4], missing2[4], missing3[4], missing4[4], missing5[4]), biobank = c(missing1[2], missing2[2], missing3[2], missing4[2], missing5[2]), barcode = c(missing1[3], missing2[3], missing3[3], missing4[3], missing5[3]), kg_number = c(missing1[1], missing2[1], missing3[1], missing4[1], missing5[1]), plate_id = rep('other', 5))
    plates = rbind(plates, df_missing, use.names=TRUE)
    table(todo$Sample %in% plates$barcode)
    # merge with plates
    todo_pheno = merge(todo, plates, by.x = 'Sample', by.y = 'barcode', all.x=T)
    # add phenotypes
    table(todo_pheno$kg_number %in% all_samples_expected$"KG-number")
    todo_pheno_diagnosis = merge(todo_pheno, all_samples_expected, by.x = 'kg_number', by.y = 'KG-number')
    # merge samples with multiple smrt cells
    todo_pheno_diagnosis_merged = data.frame()
    for (i in 1:nrow(todo_pheno_diagnosis)){
        tmp = todo_pheno_diagnosis[which(todo_pheno_diagnosis$Sample == todo_pheno_diagnosis$Sample[i]),]
        tmp = tmp[!duplicated(tmp$MovieID),]
        tmp_df = data.frame(kg_number = unique(tmp$kg_number), sample = unique(tmp$Sample), Smrt = paste(tmp$MovieID, collapse=','), yield = sum(tmp$"Total.Bases.(Gb)"), id_gwas = unique(tmp$ID_GWAS), pheno = unique(tmp$pheno))
        todo_pheno_diagnosis_merged = rbind(todo_pheno_diagnosis_merged, tmp_df)
    }
    todo_pheno_diagnosis_merged = todo_pheno_diagnosis_merged[!duplicated(todo_pheno_diagnosis_merged$sample),]
    save.image('20230316.RData')

    # ok so we need to put things together
    df_ready = merged_sequenced[, c('SMRT_CELLS', 'SAMPLE', 'COMBINED_COVERAGE', 'DIAGNOSIS', 'KG-number')]
    bbc = df_ready[is.na(df_ready$SMRT_CELLS),]; df_ready = df_ready[!is.na(df_ready$SMRT_CELLS),]
    colnames(df_ready) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
    df_single = table_occurrences_filter_info2[, c('Smrts', 'Var1', 'Cov', 'Pheno', 'KG-number')]          # modify here in case something is weird
    colnames(df_single) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
    df_nijmegen = todo_pheno_diagnosis[, c('MovieID', 'ID_GWAS', 'Total.Bases.(Gb)', 'pheno', 'kg_number')]
    colnames(df_nijmegen) = c('smrt', 'sample', 'yield', 'phenotype', 'kg_number')

    # would be nice to add the run yield for amsterdam runs
    ams_yield = data.table::fread('amsterdam_run_yield.txt', h=T, stringsAsFactors=F, dec=',')
    colnames(ams_yield) = c('smrt', 'yield')
    dim(ams_yield)
    ams_yield = ams_yield[which(ams_yield$smrt != ''),]
    ams_yield$yield = as.numeric(stringr::str_replace_all(ams_yield$yield, ',', '.'))
    # add yeild for df_ready
    nijm_logs_cp = nijm_logs
    colnames(nijm_logs_cp) = c('smrt', 'yield')
    # add a missing one
    nijm_logs_cp = rbind(nijm_logs_cp, data.frame(smrt = c('m64102e_221115_215153', 'm64102e_221113_013123', 'm64102e_221111_143527', 'm64037e_221029_002710'), yield = c(640.103, 46.344, 85.575, 46.606)))
    df_ready$yield = NA
    for (i in 1:nrow(df_ready)){
        tmp = unlist(strsplit(df_ready$smrt[i], ','))
        # search in amsterdam
        tmp_yield = ams_yield[which(ams_yield$smrt %in% tmp),]
        # search in nijmegen
        tmp_yield2 = nijm_logs_cp[which(nijm_logs_cp$smrt %in% tmp),]
        tmp_yield2 = tmp_yield2[!duplicated(tmp_yield2$smrt),]
        # combine
        tmp_combined = rbind(tmp_yield, tmp_yield2)
        # check
        if (nrow(tmp_combined) != length(tmp)){
            print(paste('Check', i))
        } else {
            df_ready$yield[i] = sum(tmp_combined$yield)
        }
    }
    # add yield for processed
    df_single$yield = NA
    for (i in 1:nrow(df_single)){
        tmp = unlist(strsplit(df_single$smrt[i], ','))
        # search in amsterdam
        tmp_yield = ams_yield[which(ams_yield$smrt %in% tmp),]
        # search in nijmegen
        tmp_yield2 = nijm_logs_cp[which(nijm_logs_cp$smrt %in% tmp),]
        tmp_yield2 = tmp_yield2[!duplicated(tmp_yield2$smrt),]
        # combine
        tmp_combined = rbind(tmp_yield, tmp_yield2)
        # check
        if (nrow(tmp_combined) != length(tmp)){
            print(paste('Check', i))
        } else {
            df_single$yield[i] = sum(tmp_combined$yield)
        }
    }

    # combine all datasets
    df_ready$type = 'ready_merged'
    df_single$type = 'to_merge_or_single'
    df_nijmegen$type = 'nijmegen_to_run'
    df_combined_with_coverage = rbind(df_ready, df_single, use.names=TRUE)

    # can we predict coverage of the nijmegen ones?
    model = lm(coverage ~ yield, data = df_combined_with_coverage)
    # fix some yield missing in df_nijmegen
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220922_051859')] = 647.275
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220624_092544')] = 579.398
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220621_142045')] = 550.855
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220816_194859.')] = 572.151
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_221002_025718')] = 34
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220625_185740')] = 458.482
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220622_234033')] = 583.096
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220627_144641')] = 527.324
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220829_124511')] = 408.805
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220830_222040')] = 440.76
    df_nijmegen$coverage = predict(model, newdata = df_nijmegen)

    # combined all together
    df_combined_all = rbind(df_combined_with_coverage, df_nijmegen)
    dim(df_combined_all)
    # we need to merge runs of the same samples again
    df_cleaned = data.frame()
    for (i in 1:nrow(df_combined_all)){
        if (!(df_combined_all$sample[i] %in% df_cleaned$sample)){
            tmp = df_combined_all[which(df_combined_all$sample == df_combined_all$sample[i]),]
            combined = data.frame(smrt = paste(tmp$smrt, collapse=','), sample = unique(tmp$sample), coverage = sum(tmp$coverage), phenotype = unique(tmp$phenotype), kg_number = unique(tmp$kg_number), yield = sum(tmp$yield), type = paste(tmp$type, collapse=','))
            df_cleaned = rbind(df_cleaned, combined)
        }
    }
    table(df_cleaned$phenotype)
    save.image('20230316.RData')
    
    # plots
    ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1")
    # low coverage
    # <12
    low_cov = df_cleaned[which(df_cleaned$coverage <12),]
    table(low_cov$phenotype)
    ggplot(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1")
    # <10
    low_cov = df_cleaned[which(df_cleaned$coverage <10),]
    dim(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),])
    table(low_cov$phenotype)
    # ad and centenarians
    ad = df_cleaned[which(df_cleaned$phenotype == 'Probable_AD'),]
    chc = df_cleaned[which(df_cleaned$phenotype == 'Centenarian'),]
    mean(ad$coverage, na.rm=T); sd(ad$coverage, na.rm=T); median(ad$coverage, na.rm=T);
    mean(chc$coverage, na.rm=T); sd(chc$coverage, na.rm=T); median(chc$coverage, na.rm=T);

    # bigger picture
    table(df_cleaned$sample %in% all_samples_expected$ID_GWAS)
    # add dates
    df_cleaned$date = NA
    for (i in 1:nrow(df_cleaned)){
        # extract dates
        tmp = paste0('20', stringr::str_split_fixed(unlist(strsplit(df_cleaned$smrt[i], ',')), '_', 3)[, 2])
        df_cleaned$date[i] = tmp[order(tmp)][1]
    }
    df_cleaned$date = as.Date(df_cleaned$date, format ='%Y%m%d')
    df_cleaned$year = NA
    df_cleaned$year[grep('2022', df_cleaned$date)] = '2022'
    df_cleaned$year[grep('2023', df_cleaned$date)] = '2023'
    df_cleaned$year[grep('2021', df_cleaned$date)] = '2021'
    year_df = data.frame(table(df_cleaned$year)); colnames(year_df) = c('Year', 'Count')
    ggplot(year_df, aes(x = Year, y = Count, fill = Year)) + geom_bar(stat = 'identity') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1")

    # write samples with low output
    low_cov = df_cleaned[which(df_cleaned$coverage <12),]
    low_cov = low_cov[which(low_cov$smrt != TRUE),]
    low_cov$type = NULL
    low_cov_smrt = unlist(strsplit(paste(low_cov$smrt, collapse = ','), ','))
    write.table(low_cov, 'samples_low_coverage.txt', quote=F, row.names=F, sep="\t")
    write.table(low_cov_smrt, 'samples_low_coverage_smrtIDs.txt', quote=F, row.names=F, sep="\t")


###############################################################
###############################################################
# Updated on 12/05
    
    # MAKE AN OVERVIEW OF THE SEQUENCING STATUS
    # LIBRARIES
    library(data.table)
    library(plyr)
    library(openxlsx)
    library(stringr)
    library(ggplot2)

    # READ ALL SAMPLES WE WERE SUPPOSED TO DO
    all_chc = fread('chc_for_pacbio_withSCORE_and_KGnumb.txt', h=T)
    all_ad = fread('20200724_final_list_ADcases_completed.txt', h=T)
    # COMBINE THEM
    colnames(all_chc) = c('ID_100plus', 'ID_GWAS', 'pheno', 'gender', 'age', 'pca_outliers_gwas', 'pbmc', 'children', 'children_partners', 'siblings', 'trios', 'brain', 'score', 'KG-number')
    colnames(all_ad) = c('I_ID', 'pheno', 'info', 'gender', 'brain_year', 'apoe', 'age_death', 'age_diag', 'last_age_available', 'brain', 'KG-number')
    # for AD cases, need to add the GWAS ID
    load('/project/holstegelab/Share/gwas_array/mapping_files/phenotypes_20211027_All.Rdata')
    tmp_pheno = pheno_final_raw[, c('I_ID', 'ID_GWAS')]
    all_ad = merge(all_ad, tmp_pheno, by = 'I_ID', all.x = T)
    # there are 3 duplicates: that is, 3 samples (AD) with 2 ID_GWAS. 2 of them make sense, as they were genotyped multiple times, 1 is more controversial as it's from the same genotyping run but different phenotypes (male/female, I_ID = 2100). Keep all duplicates in for now.
    # merge    
    all_samples_expected = rbind.fill(all_chc, all_ad)
    table(all_samples_expected$pheno)
    # we started of with 353 AD cases (3 duplicates) and 356 centenarians
    
    # READ ALL PLATES
    plates = data.frame()
    # plate 1 is in a different data, do it separately
    plate1 = read.xlsx('Plate1_QC.xlsx', sheet = 1)
    # restrict to plate, biobank, barcode and KG-number
    plate1 = plate1[, c('PlateIDtarget', 'SampleID', 'Barcode', 'Individual.ID')]
    # rename columns
    colnames(plate1) = c('plate', 'biobank', 'barcode', 'kg_number')
    # add plate id
    plate1$plate_id = 1
    # add to plates
    plates = rbind(plates, plate1)
    for (p in 2:7){
      # read sheet
      tmp = read.xlsx('All_plates_QC.xlsx', sheet = p)
      # restrict to plate, biobank, barcode and KG-number
      tmp = tmp[, c('Plate', 'Biobank.fractie', 'BC', 'KG-number')]
      # rename columns
      colnames(tmp) = c('plate', 'biobank', 'barcode', 'kg_number')
      # add plate id
      tmp$plate_id = p
      # add to plates
      plates = rbind(plates, tmp)
    }
    table(plates$plate_id)
    # the plates contain 653 samples in total
    
    # NOW READ ALL SEQUENCED SAMPLES BASED ON THE COVERAGE SUMMARY FILE
    merged_sequenced = fread('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', h=T, stringsAsFactors=F, sep="\t", dec = ',')
    merged_sequenced = merged_sequenced[!is.na(merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[!duplicated(merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[order(-merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[!duplicated(merged_sequenced$SAMPLE),]
    # also read blood-brain-child samples
    bbc_samples = fread('/project/holstegelab/Software/snakemake_pipeline/sample_check_data/20210615_overview_genetics.txt', h=T, stringsAsFactors=F)
    bbc_samples = bbc_samples[!is.na(bbc_samples$ID_PACBIO),]
    # modify here
    bbc_samples = bbc_samples[, c('ID_GWAS')]
    bbc_samples$DIAGNOSIS = 'Blood_Brain_Child'
    colnames(bbc_samples) = c('SAMPLE', 'DIAGNOSIS')
    merged_sequenced = rbind.fill(merged_sequenced, bbc_samples)
    table(merged_sequenced$DIAGNOSIS)
    # in total, ready for analysis there are 69 AD cases + 84 Centenarians + 10 more centenarians (blood-brain-child project) + 10 children (blood-brain-child project)
    
    # NOW LABEL THE SAMPLES SEQUENCED IN ALL_SAMPLES WE WERE SUPPOSED TO SEQUENCE
    table(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)
    # 163/173 are OK -- let's check the remaining 10
    non_matching = merged_sequenced[which(!(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)),]
    # add phenotypes
    non_matching_pheno = merge(non_matching, pheno_final_raw, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    # it's ok -- these 10 are the children of the centenarians
    # So all samples we sequenced until now were in the list of samples to be sequenced
    
    # NOW CHECK WITH RESPECT TO THE PLATES
    # add kg-number to the sequenced samples
    kg_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
    merged_sequenced = merge(merged_sequenced, kg_info, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    table(merged_sequenced$DIAGNOSIS)
    # merge sequenced samples with plate information
    merged_sequenced_plate = merge(merged_sequenced, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    # add phenotypes
    tmp_pheno = pheno_final_raw[, c('ID_GWAS', 'diagnosis', 'ID_100plus', 'I_ID')]
    merged_sequenced_plate = merge(merged_sequenced_plate, tmp_pheno, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    table(merged_sequenced_plate$DIAGNOSIS)
    # there are duplicates, we should clean them up
    dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
    merged_sequenced_plate = merged_sequenced_plate[order(-merged_sequenced_plate$COMBINED_COVERAGE),]
    merged_sequenced_plate = merged_sequenced_plate[!duplicated(merged_sequenced_plate$SAMPLE),] 
    dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
    # find samples wiith missing plate id
    merged_sequenced_plate[is.na(merged_sequenced_plate$plate_id),]
    table(merged_sequenced_plate$plate_id, exclude=F)
    merged_sequenced_plate$SMRT_CELL_N = factor(merged_sequenced_plate$SMRT_CELL_N)
    # 15 in total: 10 are centenarians from the blood-brain-child project, 4 are additional centenarians (not sure why they are not in the list) + 1 AD (ND-001200)
    
    # PLOT THE COMBINED_COVERAGE OF SAMPLES
    pdf('coverage_merged_hq_samples.pdf', height = 7, width = 10)
    ggplot(merged_sequenced_plate[which(merged_sequenced_plate$DIAGNOSIS != 'Blood_Brain_Child'),], aes(x = DIAGNOSIS, y = COMBINED_COVERAGE, fill = DIAGNOSIS)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = SMRT_CELL_N), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Combined coverage of HiFi reads')
    dev.off()
    
    # LOOK AT SMRT CELLS >600GB YIELD
    smrt_coverage_p1 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-02-17_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    smrt_coverage_p2 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-03-10_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    smrt_coverage_p3 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-05-12_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    # subset columns
    smrt_coverage_p3 = smrt_coverage_p3[, c('movie_id', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    smrt_coverage_p2 = smrt_coverage_p2[, c('movie_id', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    smrt_coverage_p1 = smrt_coverage_p1[, c('SMRT_ID', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    colnames(smrt_coverage_p2)[1] = 'SMRT_ID'
    colnames(smrt_coverage_p3)[1] = 'SMRT_ID'
    # combine
    smrt_coverage_all = rbind(smrt_coverage_p1, smrt_coverage_p2, smrt_coverage_p3)
    # remove duplicates based on smrt id
    dups = smrt_coverage_all[duplicated(smrt_coverage_all$SMRT_ID),]
    smrt_coverage_all = smrt_coverage_all[!duplicated(smrt_coverage_all$SMRT_ID),]
    dim(smrt_coverage_all)
    # 706 good smrt cells sequenced thus far
    # ok now exclude samples combined (from merged_sequenced_plate)
    smrt_coverage_all_notcombined = smrt_coverage_all[which(!(smrt_coverage_all$ID_GWAS %in% merged_sequenced_plate$SAMPLE)),]
    # also exclude samples from blood-brain-child project
    smrt_coverage_all_notcombined = smrt_coverage_all_notcombined[which(!(smrt_coverage_all_notcombined$ID_GWAS %in% bbc_samples$SAMPLE)),]
    table_occurrences = data.frame(table(smrt_coverage_all_notcombined$ID_GWAS))
    table_occurrences = table_occurrences[order(-table_occurrences$Freq),]
    head(table_occurrences)
    # calculate total coverage of these
    table_occurrences$Cov = NA
    table_occurrences$Pheno = NA
    table_occurrences$Smrts = NA
    for (i in 1:nrow(table_occurrences)){
      tmp = smrt_coverage_all_notcombined[which(smrt_coverage_all_notcombined$ID_GWAS == table_occurrences$Var1[i]),]
      table_occurrences$Cov[i] = sum(tmp$GLOBAL_COVERAGE)
      table_occurrences$Pheno[i] = unique(tmp$diagnosis)
      table_occurrences$Smrts[i] = paste(tmp$SMRT_ID, collapse = ',')
    }
    table_occurrences = table_occurrences[order(-table_occurrences$Cov),]
    head(table_occurrences, 10)

    # get most updated list of samples merged
    samples_merged_chc = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/centenarian/*merged.hifi.chm13.bam', intern=T)
    samples_merged_ad = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/ad/*merged.hifi.chm13.bam', intern=T)
    df_chc = data.frame(sample = samples_merged_chc, pheno = 'chc')
    df_ad = data.frame(sample = samples_merged_ad, pheno = 'ad')
    df_combined = rbind(df_chc, df_ad)
    for (i in 1:nrow(df_combined)){
      tmp = str_replace_all(strsplit(df_combined$sample[i], '/')[[1]][length(strsplit(df_combined$sample[i], '/')[[1]])], '.merged.hifi.chm13.bam', '')
      df_combined$sample[i] = tmp
    }
    # exclude from the table of non combined any of these actually combined
    dim(table_occurrences)
    table_occurrences_filter = table_occurrences[which(!(table_occurrences$Var1 %in% df_combined$sample)), ]
    dim(table_occurrences_filter)

    # then take out the samples with coverage >12
    coverage_enouth = table_occurrences_filter[which(table_occurrences_filter$Cov >= 12), ]
    dim(coverage_enouth)
    coverage_enouth$Freq = factor(coverage_enouth$Freq)
    coverage_enouth$Pheno[which(coverage_enouth$Pheno == 'Control_100plus')] = 'Centenarian'
    
    # plot all samples with actual data
    pdf('all_samples_with_data.pdf', height = 7, width = 10)
    table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno == 'Control_100plus')] = 'Centenarian'
    table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno == 'AD_othersite')] = 'Probable_AD'
    table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno %in% c('NIID', 'Dementie_anders'))] = 'Anke'
    table_occurrences_filter = table_occurrences_filter[which(table_occurrences_filter$Pheno != 'family_100plus'),]
    table_occurrences_filter$Freq = factor(table_occurrences_filter$Freq)
    ggplot(table_occurrences_filter, aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + ylim(0, 45) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    dev.off()
    # plot samples with enough coverage
    pdf('samples_with_enough_coverage_data.pdf', height = 7, width = 10)
    ggplot(coverage_enouth, aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + ylim(10, 45) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    dev.off()
    # plot samples without enough coverage
    pdf('samples_without_enough_coverage_data.pdf', height = 7, width = 10)
    coverage_low = table_occurrences_filter[which(table_occurrences_filter$Cov <12),]
    coverage_low$Freq = factor(coverage_low$Freq)
    ggplot(coverage_low, aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    dev.off()
    # same but with histogram
    ggplot(coverage_low, aes(x = Cov)) + geom_histogram()

    # also overlap with plates
    plates_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
    coverage_enouth = merge(coverage_enouth, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
    coverage_enouth_plate = merge(coverage_enouth, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    # check plates
    table(merged_sequenced_plate$plate_id, exclude=F)
    table(coverage_enouth_plate$plate_id, exclude=F)
    # all samples with data
    table_occurrences_filter_info = merge(table_occurrences_filter, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
    table_occurrences_filter_info_plate = merge(table_occurrences_filter_info, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    table(table_occurrences_filter_info_plate$plate_id, exclude=F)

    # next go to the Nijmegen stats
    nijm_logs = fread('Nijmegen_Log.txt', h=T, stringsAsFactors=F, dec=',')
    # exclude runa we have processed
    nijm_logs_missing = nijm_logs[which(!(nijm_logs$MovieID %in% smrt_coverage_all$SMRT_ID)),]
    # n=22 (makes sense, 1 run from june, some from October and November -- probably they still have to come)
    # plot yield and date
    nijm_logs_missing$date = str_split_fixed(nijm_logs_missing$MovieID, '_', 3)[, 2]
    nijm_logs_missing$date = paste0('20', nijm_logs_missing$date)
    nijm_logs_missing$date = as.Date(nijm_logs_missing$date, "%Y%m%d")
    colnames(nijm_logs_missing) = c('Smrt', 'Yield', 'Date')
    ggplot(nijm_logs_missing, aes(x=Date, y=Yield, color=Yield)) + geom_point(stat = 'identity', size=3) + ggtitle('Missing runs from Nijmegen only based on Log files') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
    # see how many runs with yield >600gb we have
    dim(nijm_logs_missing[which(nijm_logs_missing$Yield >= 600),])
    # maybe plot Yield over time (all smrt cells ever?)
    all_run_logs = fread('All_runs_Logs.txt', h=T, stringsAsFactors=F, dec=',')
    all_run_logs$date = str_split_fixed(all_run_logs$MovieID, '_', 3)[, 2]
    all_run_logs$date = paste0('20', all_run_logs$date)
    all_run_logs$date = as.Date(all_run_logs$date, "%Y%m%d")
    colnames(all_run_logs) = c('Smrt', 'Yield', 'Site', 'Date')
    all_run_logs = all_run_logs[!is.na(all_run_logs$Yield),]
    all_run_logs = all_run_logs[order(all_run_logs$Date),]
    ggplot(all_run_logs, aes(x=Date, y=Yield, color=Site)) + geom_point(stat = 'identity', size=3, alpha = 0.5) + geom_smooth(method = 'loess') + ggtitle('SMRT Yield over time and sequencing centers') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
    
    #################################### -- Final integrating all things
    # AIM IS TO GENERATE AN OVERVIEW AT THE SAMPLE LEVEL
    # This is based on the merged samples
    merged_sequenced
    # This is based on the single smrt cells
    table_occurrences_filter
    # check whether we expected to sequence these
    tofix = table_occurrences_filter[which(!(table_occurrences_filter$Var1 %in% all_samples_expected$ID_GWAS)),]
    # 2 centenarians were (apparently) not supposed to be sequenced (nor sure why, nothing weird. maybe sample swap??), and 1 AD (EMC_NEURO study, probably also Anke/Annemieke). 3 Anke's samples are OK.
    table_occurrences_filter_info = merge(table_occurrences_filter, pheno_final_raw, by.x = 'Var1', by.y = 'ID_GWAS')
    table_occurrences_filter_info2 = merge(table_occurrences_filter_info, all_samples_expected, by.x = 'Var1', by.y = 'ID_GWAS', all.x=T)
    
    # then we need to go to Nijmegen data
    new_nijm = read.xlsx('Nijmegen_Sample_ID_and_RunIDs_Overview.xlsx')
    # exclude samples from the blood-brain-child project
    new_nijm = new_nijm[which(new_nijm$Project != 'BBC'),]
    # from this new list, exclude smrt cells already processed (in table occurrences)
    new_nijm$Done = NA
    for (i in 1:nrow(new_nijm)){
      # get run id
      tmp_smrt = new_nijm$MovieID[i]
      # check if was merged
      was_merged = merged_sequenced[grep(tmp_smrt, merged_sequenced$SMRT_CELLS),]
      was_proces = table_occurrences_filter[grep(tmp_smrt, table_occurrences_filter$Smrts),]
      if (nrow(was_merged) >0){
        new_nijm$Done[i] = 'merged'
      } else if (nrow(was_proces) >0){
        new_nijm$Done[i] = 'proces'
      } else{
        new_nijm$Done[i] = 'todo'
      }
    }
    table(new_nijm$Done, exclude=F)
    # then take only todo
    todo = new_nijm[which(new_nijm$Done == 'todo'),]
    todo$Sample = stringr::str_replace_all(todo$Sample, '-02', '')
    todo$Sample = stringr::str_replace_all(todo$Sample, '_19kb_30-07-2021', '')
    
    # read run stats
    run_stats = read.xlsx('Copy_nijmegen_runs_from_okt.xlsx', sheet = 2)
    table(todo$MovieID %in% run_stats$context)
    not_matching = todo[which(!todo$MovieID %in% run_stats$context),]
    todo = merge(todo, run_stats, by.x = 'MovieID', by.y = 'context', all.x=T)
    # 16 do not match between the run statistics from Nijmegen: this can be because the run statistics (run_stats) are only from october onwards, while the other file (todo) also includes other months, so that's ok
    table(todo$Sample %in% plates$barcode)
    # 7 samples sequenced were not in the plates -- maybe some change we did afterwards?
    mismatch = todo[which(!(todo$Sample %in% plates$barcode)),]
    # solve these 7
    missing1 = c('KG-011524', '17R1457', 'FR34137283', '202204011P1V', 'D01')
    missing2 = c('ND-001200', '11D01145', 'FR34137218', '202204011P1V', 'A01')
    missing3 = c('ND-000882', '10D04356', 'FR34137285', '202204011P1V', 'G04')
    missing4 = c('ND-000025', '09D04551', 'FR34137227', '202204011P1V', 'B06')
    missing5 = c('KG-011559', '17R1369', 'FR34137302', '202204011P1V', 'H01')
    df_missing = data.frame(plate = c(missing1[4], missing2[4], missing3[4], missing4[4], missing5[4]), biobank = c(missing1[2], missing2[2], missing3[2], missing4[2], missing5[2]), barcode = c(missing1[3], missing2[3], missing3[3], missing4[3], missing5[3]), kg_number = c(missing1[1], missing2[1], missing3[1], missing4[1], missing5[1]), plate_id = rep('other', 5))
    plates = rbind(plates, df_missing, use.names=TRUE)
    table(todo$Sample %in% plates$barcode)
    # merge with plates
    todo_pheno = merge(todo, plates, by.x = 'Sample', by.y = 'barcode', all.x=T)
    # add phenotypes
    table(todo_pheno$kg_number %in% all_samples_expected$"KG-number")
    todo_pheno_diagnosis = merge(todo_pheno, all_samples_expected, by.x = 'kg_number', by.y = 'KG-number')
    # merge samples with multiple smrt cells
    todo_pheno_diagnosis_merged = data.frame()
    for (i in 1:nrow(todo_pheno_diagnosis)){
      tmp = todo_pheno_diagnosis[which(todo_pheno_diagnosis$Sample == todo_pheno_diagnosis$Sample[i]),]
      tmp = tmp[!duplicated(tmp$MovieID),]
      tmp_df = data.frame(kg_number = unique(tmp$kg_number), sample = unique(tmp$Sample), Smrt = paste(tmp$MovieID, collapse=','), yield = sum(tmp$"Total.Bases.(Gb)"), id_gwas = unique(tmp$ID_GWAS), pheno = unique(tmp$pheno))
      todo_pheno_diagnosis_merged = rbind(todo_pheno_diagnosis_merged, tmp_df)
    }
    todo_pheno_diagnosis_merged = todo_pheno_diagnosis_merged[!duplicated(todo_pheno_diagnosis_merged$sample),]
    
    ######################################## 
    # SO NOW WE NEED TO PUT THINGS TOGETHER
    df_ready = merged_sequenced[, c('SMRT_CELLS', 'SAMPLE', 'COMBINED_COVERAGE', 'DIAGNOSIS', 'KG-number')]
    bbc = df_ready[is.na(df_ready$SMRT_CELLS),]; df_ready = df_ready[!is.na(df_ready$SMRT_CELLS),]
    colnames(df_ready) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
    df_single = table_occurrences_filter_info2[, c('Smrts', 'Var1', 'Cov', 'Pheno', 'KG-number')]          # modify here in case something is weird
    colnames(df_single) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
    df_nijmegen = todo_pheno_diagnosis[, c('MovieID', 'ID_GWAS', 'Total.Bases.(Gb)', 'pheno', 'kg_number')]
    colnames(df_nijmegen) = c('smrt', 'sample', 'yield', 'phenotype', 'kg_number')

    # would be nice to add the run yield for amsterdam runs
    ams_yield = data.table::fread('amsterdam_run_yield_20230512.txt', h=T, stringsAsFactors=F, dec=',')
    colnames(ams_yield) = c('smrt', 'yield')
    dim(ams_yield)
    ams_yield = ams_yield[which(ams_yield$smrt != ''),]
    ams_yield$yield = as.numeric(stringr::str_replace_all(ams_yield$yield, ',', '.'))
    
    # nijmegen needs to be combined
    nijm_logs_cp = nijm_logs
    colnames(nijm_logs_cp) = c('smrt', 'yield')
    run_stat_nj = run_stats[, c('context', 'Total.Bases.(Gb)')]
    run_stat_nj = run_stat_nj[which(!run_stat_nj$context %in% nijm_logs_cp$smrt),]
    colnames(run_stat_nj) = c('smrt', 'yield')
    nijm_logs_cp = rbind(nijm_logs_cp, run_stat_nj)

    # add yield for df_ready     
    # add a missing one
    nijm_logs_cp = rbind(nijm_logs_cp, data.frame(smrt = c('m64102e_221115_215153', 'm64102e_221113_013123', 'm64102e_221111_143527', 'm64037e_221029_002710'), yield = c(640.103, 46.344, 85.575, 46.606)))
    df_ready$yield = NA
    for (i in 1:nrow(df_ready)){
      tmp = unlist(strsplit(df_ready$smrt[i], ','))
      # search in amsterdam
      tmp_yield = ams_yield[which(ams_yield$smrt %in% tmp),]
      # search in nijmegen
      tmp_yield2 = nijm_logs_cp[which(nijm_logs_cp$smrt %in% tmp),]
      tmp_yield2 = tmp_yield2[!duplicated(tmp_yield2$smrt),]
      # combine
      tmp_combined = rbind(tmp_yield, tmp_yield2)
      # check
      if (nrow(tmp_combined) != length(tmp)){
        print(paste('Check', i))
      } else {
        df_ready$yield[i] = sum(tmp_combined$yield)
      }
    }
    # add yield for processed
    df_single$yield = NA
    for (i in 1:nrow(df_single)){
      tmp = unlist(strsplit(df_single$smrt[i], ','))
      # search in amsterdam
      tmp_yield = ams_yield[which(ams_yield$smrt %in% tmp),]
      # search in nijmegen
      tmp_yield2 = nijm_logs_cp[which(nijm_logs_cp$smrt %in% tmp),]
      tmp_yield2 = tmp_yield2[!duplicated(tmp_yield2$smrt),]
      # combine
      tmp_combined = rbind(tmp_yield, tmp_yield2)
      # check
      if (nrow(tmp_combined) != length(tmp)){
        print(paste('Check', i))
      } else {
        df_single$yield[i] = sum(tmp_combined$yield)
      }
    }
    
    # combine all datasets
    df_ready$type = 'ready_merged'
    df_single$type = 'to_merge_or_single'
    df_nijmegen$type = 'nijmegen_to_run'
    df_combined_with_coverage = rbind(df_ready, df_single, use.names=TRUE)
    
    # show correlation between yield and coverage
    ggplot(data = df_combined_with_coverage, aes(x = coverage, y = yield)) + geom_point(stat = 'identity', size = 3) + geom_smooth() + xlab('Coverage of HiFi data') + ylab('Total Yield (GB)')

    # can we predict coverage of the nijmegen ones?
    model = lm(coverage ~ yield, data = df_combined_with_coverage)
    # fix some yield missing in df_nijmegen
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220922_051859')] = 647.275
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220624_092544')] = 579.398
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220621_142045')] = 550.855
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220816_194859.')] = 572.151
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_221002_025718')] = 34
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220625_185740')] = 458.482
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220622_234033')] = 583.096
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220627_144641')] = 527.324
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220829_124511')] = 408.805
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220830_222040')] = 440.76
    df_nijmegen$coverage = predict(model, newdata = df_nijmegen)
    
    # combined all together
    df_combined_all = rbind(df_combined_with_coverage, df_nijmegen)
    dim(df_combined_all)
    # we need to merge runs of the same samples again
    df_cleaned = data.frame()
    for (i in 1:nrow(df_combined_all)){
      if (!(df_combined_all$sample[i] %in% df_cleaned$sample)){
        tmp = df_combined_all[which(df_combined_all$sample == df_combined_all$sample[i]),]
        combined = data.frame(smrt = paste(tmp$smrt, collapse=','), sample = unique(tmp$sample), coverage = sum(tmp$coverage, na.rm=T), phenotype = unique(tmp$phenotype), kg_number = unique(tmp$kg_number), yield = sum(tmp$yield, na.rm=T), type = paste(tmp$type, collapse=','))
        df_cleaned = rbind(df_cleaned, combined)
      }
    }
    table(df_cleaned$phenotype)
    df_cleaned = df_cleaned[which(df_cleaned$smrt != 'TRUE'),]
    table(df_cleaned$phenotype)
    
    # plots
    # all samples, Amsterdam and Nijmegen together
    pdf('all_samples_with_prediction.pdf', height=7, width=10)
    ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
    dev.off()

    # same plot but color good and bad samples
    df_cleaned$Status = ifelse(df_cleaned$coverage >=12, 'Good', 'Re-do')
    pdf('all_samples_with_prediction_status.pdf', height=7, width=10)
    ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2, aes(col = Status)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
    dev.off()

    # low coverage
    # <12
    pdf('redo_with_predicted.pdf', height = 7, width = 10)
    low_cov = df_cleaned[which(df_cleaned$coverage <12),]
    table(low_cov$phenotype)
    ggplot(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data') + 
    dev.off()
    # also histogram
    pdf('redo_with_predicted_hist.pdf', height = 7, width = 10)
    ggplot(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = coverage, fill = phenotype)) + geom_histogram() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Coverage of HiFi data') + ylab('Samples count')
    dev.off()
    save.image('20230512.RData')

    # output the list of samples with low coverage
    low_cov_ad_chc = low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),]
    table(low_cov_ad_chc$kg_number %in% plates$kg_number)
    low_cov_ad_chc_info = merge(low_cov_ad_chc, plates, by = 'kg_number', all.x = T)
    low_cov_ad_chc_info = low_cov_ad_chc_info[!duplicated(low_cov_ad_chc_info$kg_number),]
    write.table(low_cov_ad_chc_info, '20230515_low_coverage_samples.txt', quote=F, row.names=F, sep="\t", dec=',')
    save.image('20230512.RData')

###############################################################
###############################################################
# Updated on 16/06
    # MAKE AN OVERVIEW OF THE SEQUENCING STATUS
    # LIBRARIES
    library(data.table)
    library(plyr)
    library(openxlsx)
    library(stringr)
    library(ggplot2)

    # READ ALL SAMPLES WE WERE SUPPOSED TO DO
    all_chc = fread('chc_for_pacbio_withSCORE_and_KGnumb.txt', h=T)
    all_ad = fread('20200724_final_list_ADcases_completed.txt', h=T)
    # COMBINE THEM
    colnames(all_chc) = c('ID_100plus', 'ID_GWAS', 'pheno', 'gender', 'age', 'pca_outliers_gwas', 'pbmc', 'children', 'children_partners', 'siblings', 'trios', 'brain', 'score', 'KG-number')
    colnames(all_ad) = c('I_ID', 'pheno', 'info', 'gender', 'brain_year', 'apoe', 'age_death', 'age_diag', 'last_age_available', 'brain', 'KG-number')
    # for AD cases, need to add the GWAS ID
    load('/project/holstegelab/Share/gwas_array/mapping_files/phenotypes_20211027_All.Rdata')
    tmp_pheno = pheno_final_raw[, c('I_ID', 'ID_GWAS')]
    all_ad = merge(all_ad, tmp_pheno, by = 'I_ID', all.x = T)
    # there are 3 duplicates: that is, 3 samples (AD) with 2 ID_GWAS. 2 of them make sense, as they were genotyped multiple times, 1 is more controversial as it's from the same genotyping run but different phenotypes (male/female, I_ID = 2100). Keep all duplicates in for now.
    # merge
    all_samples_expected = rbind.fill(all_chc, all_ad)
    table(all_samples_expected$pheno)
    # we started of with 353 AD cases (3 duplicates) and 356 centenarians
    
    # READ ALL PLATES
    plates = data.frame()
    # plate 1 is in a different data, do it separately
    plate1 = read.xlsx('Plate1_QC.xlsx', sheet = 1)
    # restrict to plate, biobank, barcode and KG-number
    plate1 = plate1[, c('PlateIDtarget', 'SampleID', 'Barcode', 'Individual.ID')]
    # rename columns
    colnames(plate1) = c('plate', 'biobank', 'barcode', 'kg_number')
    # add plate id
    plate1$plate_id = 1
    # add to plates
    plates = rbind(plates, plate1)
    for (p in 2:7){
      # read sheet
      tmp = read.xlsx('All_plates_QC.xlsx', sheet = p)
      # restrict to plate, biobank, barcode and KG-number
      tmp = tmp[, c('Plate', 'Biobank.fractie', 'BC', 'KG-number')]
      # rename columns
      colnames(tmp) = c('plate', 'biobank', 'barcode', 'kg_number')
      # add plate id
      tmp$plate_id = p
      # add to plates
      plates = rbind(plates, tmp)
    }
    table(plates$plate_id)
    # the plates contain 653 samples in total
    
    # NOW READ ALL SEQUENCED SAMPLES BASED ON THE COVERAGE SUMMARY FILE
    merged_sequenced = fread('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', h=T, stringsAsFactors=F, sep="\t", dec = ',')
    merged_sequenced = merged_sequenced[!is.na(merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[order(-merged_sequenced$COMBINED_COVERAGE),]
    merged_sequenced = merged_sequenced[!duplicated(merged_sequenced$SAMPLE),]
    # also read blood-brain-child samples
    bbc_samples = fread('/project/holstegelab/Software/snakemake_pipeline/sample_check_data/20210615_overview_genetics.txt', h=T, stringsAsFactors=F)
    bbc_samples = bbc_samples[!is.na(bbc_samples$ID_PACBIO),]
    # modify here
    bbc_samples = bbc_samples[, c('ID_GWAS')]
    bbc_samples$DIAGNOSIS = 'Blood_Brain_Child'
    colnames(bbc_samples) = c('SAMPLE', 'DIAGNOSIS')
    merged_sequenced = rbind.fill(merged_sequenced, bbc_samples)
    table(merged_sequenced$DIAGNOSIS)
    # in total, ready for analysis there are 81 AD cases + 104 Centenarians + 10 more centenarians (blood-brain-child project) + 10 children (blood-brain-child project)
    
    # NOW LABEL THE SAMPLES SEQUENCED IN ALL_SAMPLES WE WERE SUPPOSED TO SEQUENCE
    table(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)
    # 179/189 are OK -- let's check the remaining 10
    non_matching = merged_sequenced[which(!(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)),]
    # add phenotypes
    non_matching_pheno = merge(non_matching, pheno_final_raw, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    # it's ok -- these 10 are the children of the centenarians
    # So all samples we sequenced until now were in the list of samples to be sequenced
    
    # NOW CHECK WITH RESPECT TO THE PLATES
    # add kg-number to the sequenced samples
    kg_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
    merged_sequenced = merge(merged_sequenced, kg_info, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    table(merged_sequenced$DIAGNOSIS)
    # merge sequenced samples with plate information
    merged_sequenced_plate = merge(merged_sequenced, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    # add phenotypes
    tmp_pheno = pheno_final_raw[, c('ID_GWAS', 'diagnosis', 'ID_100plus', 'I_ID')]
    merged_sequenced_plate = merge(merged_sequenced_plate, tmp_pheno, by.x = 'SAMPLE', by.y = 'ID_GWAS')
    table(merged_sequenced_plate$DIAGNOSIS)
    # there are duplicates, we should clean them up
    dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
    merged_sequenced_plate = merged_sequenced_plate[order(-merged_sequenced_plate$COMBINED_COVERAGE),]
    merged_sequenced_plate = merged_sequenced_plate[!duplicated(merged_sequenced_plate$SAMPLE),] 
    dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
    # find samples wiith missing plate id
    merged_sequenced_plate[is.na(merged_sequenced_plate$plate_id),]
    table(merged_sequenced_plate$plate_id, exclude=F)
    merged_sequenced_plate$SMRT_CELL_N = factor(merged_sequenced_plate$SMRT_CELL_N)
    # 15 in total: 10 are centenarians from the blood-brain-child project, 4 are additional centenarians (not sure why they are not in the list) + 1 AD (ND-001200)
    
    # PLOT THE COMBINED_COVERAGE OF SAMPLES
    pdf('coverage_merged_hq_samples.pdf', height = 7, width = 10)
    ggplot(merged_sequenced_plate[which(merged_sequenced_plate$DIAGNOSIS != 'Blood_Brain_Child'),], aes(x = DIAGNOSIS, y = COMBINED_COVERAGE, fill = DIAGNOSIS)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(shape = SMRT_CELL_N), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Combined coverage of HiFi reads')
    dev.off()
    
    # LOOK AT SMRT CELLS >600GB YIELD
    smrt_coverage_p1 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-02-17_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    smrt_coverage_p2 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-03-10_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    smrt_coverage_p3 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-05-12_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    smrt_coverage_p3 = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-06-16_freeze_sequencing_stats.txt', h=T, stringsAsFactors=F, dec = ',')
    # subset columns
    smrt_coverage_p3 = smrt_coverage_p3[, c('movie_id', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    smrt_coverage_p2 = smrt_coverage_p2[, c('movie_id', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    smrt_coverage_p1 = smrt_coverage_p1[, c('SMRT_ID', 'ID_GWAS', 'original_ID', 'PERC_HOMOLOGY', 'ID_100plus', 'GLOBAL_COVERAGE', 'diagnosis')]
    colnames(smrt_coverage_p2)[1] = 'SMRT_ID'
    colnames(smrt_coverage_p3)[1] = 'SMRT_ID'
    # combine
    smrt_coverage_all = rbind(smrt_coverage_p1, smrt_coverage_p2, smrt_coverage_p3)
    # remove duplicates based on smrt id
    dups = smrt_coverage_all[duplicated(smrt_coverage_all$SMRT_ID),]
    smrt_coverage_all = smrt_coverage_all[!duplicated(smrt_coverage_all$SMRT_ID),]
    dim(smrt_coverage_all)
    # 720 good smrt cells sequenced thus far --> make output to update my Excel
    write.table(smrt_coverage_all, '/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-06-16_combined_coverage_all.txt', quote=F, row.names=F, sep='\t', dec=',')
    # ok now exclude samples combined (from merged_sequenced_plate)
    smrt_coverage_all_notcombined = smrt_coverage_all[which(!(smrt_coverage_all$ID_GWAS %in% merged_sequenced_plate$SAMPLE)),]
    # also exclude samples from blood-brain-child project
    smrt_coverage_all_notcombined = smrt_coverage_all_notcombined[which(!(smrt_coverage_all_notcombined$ID_GWAS %in% bbc_samples$SAMPLE)),]
    table_occurrences = data.frame(table(smrt_coverage_all_notcombined$ID_GWAS))
    table_occurrences = table_occurrences[order(-table_occurrences$Freq),]
    head(table_occurrences)
    # calculate total coverage of these
    table_occurrences$Cov = NA
    table_occurrences$Pheno = NA
    table_occurrences$Smrts = NA
    for (i in 1:nrow(table_occurrences)){
      tmp = smrt_coverage_all_notcombined[which(smrt_coverage_all_notcombined$ID_GWAS == table_occurrences$Var1[i]),]
      table_occurrences$Cov[i] = sum(tmp$GLOBAL_COVERAGE)
      table_occurrences$Pheno[i] = unique(tmp$diagnosis)
      table_occurrences$Smrts[i] = paste(tmp$SMRT_ID, collapse = ',')
    }
    table_occurrences = table_occurrences[order(-table_occurrences$Cov),]
    head(table_occurrences, 10)

    # get most updated list of samples merged
    samples_merged_chc = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/centenarian/*merged.hifi.chm13.bam', intern=T)
    samples_merged_ad = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/ad/*merged.hifi.chm13.bam', intern=T)
    df_chc = data.frame(sample = samples_merged_chc, pheno = 'chc')
    df_ad = data.frame(sample = samples_merged_ad, pheno = 'ad')
    df_combined = rbind(df_chc, df_ad)
    for (i in 1:nrow(df_combined)){
      tmp = str_replace_all(strsplit(df_combined$sample[i], '/')[[1]][length(strsplit(df_combined$sample[i], '/')[[1]])], '.merged.hifi.chm13.bam', '')
      df_combined$sample[i] = tmp
    }
    # exclude from the table of non combined any of these actually combined
    dim(table_occurrences)
    table_occurrences_filter = table_occurrences[which(!(table_occurrences$Var1 %in% df_combined$sample)), ]
    dim(table_occurrences_filter)
    # there is 1 centenrian duplicated, the rest is ok
    table_occurrences_filter = table_occurrences_filter[which(table_occurrences_filter$Var1 != '3196128'),]

    # then take out the samples with coverage >12
    coverage_enouth = table_occurrences_filter[which(table_occurrences_filter$Cov >= 12), ]
    dim(coverage_enouth)
    coverage_enouth$Freq = factor(coverage_enouth$Freq)
    coverage_enouth$Pheno[which(coverage_enouth$Pheno == 'Control_100plus')] = 'Centenarian'
    coverage_enouth$Pheno[which(coverage_enouth$Pheno == 'MCI')] = 'Probable_AD'
    
    # plot all samples with actual data
    pdf('all_samples_with_data.pdf', height = 7, width = 10)
    table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno == 'Control_100plus')] = 'Centenarian'
    table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno %in% c('AD_othersite', 'MCI'))] = 'Probable_AD'
    table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno %in% c('NIID', 'Dementie_anders'))] = 'Anke'
    table_occurrences_filter = table_occurrences_filter[which(table_occurrences_filter$Pheno != 'family_100plus'),]
    table_occurrences_filter$Freq = factor(table_occurrences_filter$Freq)
    ggplot(table_occurrences_filter, aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + ylim(0, 45) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    dev.off()
    # plot samples with enough coverage
    pdf('samples_with_enough_coverage_data.pdf', height = 7, width = 10)
    ggplot(coverage_enouth, aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + ylim(10, 45) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    dev.off()
    # plot samples without enough coverage
    pdf('samples_without_enough_coverage_data.pdf', height = 7, width = 10)
    coverage_low = table_occurrences_filter[which(table_occurrences_filter$Cov <12),]
    coverage_low$Freq = factor(coverage_low$Freq)
    ggplot(coverage_low, aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
    dev.off()

    # also overlap with plates
    plates_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
    coverage_enouth = merge(coverage_enouth, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
    coverage_enouth_plate = merge(coverage_enouth, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    # check plates
    table(merged_sequenced_plate$plate_id, exclude=F)
    table(coverage_enouth_plate$plate_id, exclude=F)
    # all samples with data
    table_occurrences_filter_info = merge(table_occurrences_filter, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
    table_occurrences_filter_info_plate = merge(table_occurrences_filter_info, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
    table(table_occurrences_filter_info_plate$plate_id, exclude=F)

    # next go to the Nijmegen stats
    nijm_logs = fread('Nijmegen_Log.txt', h=T, stringsAsFactors=F, dec=',')
    # exclude runa we have processed
    nijm_logs_missing = nijm_logs[which(!(nijm_logs$MovieID %in% smrt_coverage_all$SMRT_ID)),]
    # n=22 (makes sense, 1 run from june, some from October and November -- probably they still have to come)
    # plot yield and date
    nijm_logs_missing$date = str_split_fixed(nijm_logs_missing$MovieID, '_', 3)[, 2]
    nijm_logs_missing$date = paste0('20', nijm_logs_missing$date)
    nijm_logs_missing$date = as.Date(nijm_logs_missing$date, "%Y%m%d")
    colnames(nijm_logs_missing) = c('Smrt', 'Yield', 'Date')
    # add the very new data
    nijm_logs_newer = fread('Nijmegen_Log_new.txt', h=T, stringsAsFactors=F, dec=',')
    # exclude runa we have processed
    nijm_logs_newer_missing = nijm_logs_newer[which(!(nijm_logs_newer$MovieID %in% smrt_coverage_all$SMRT_ID)),]
    # n=22 (makes sense, 1 run from june, some from October and November -- probably they still have to come)
    # plot yield and date
    nijm_logs_newer_missing$date = str_split_fixed(nijm_logs_newer_missing$MovieID, '_', 3)[, 2]
    nijm_logs_newer_missing$date = paste0('20', nijm_logs_newer_missing$date)
    nijm_logs_newer_missing$date = as.Date(nijm_logs_newer_missing$date, "%Y%m%d")
    colnames(nijm_logs_newer_missing) = c('Smrt', 'Yield', 'Date')
    # combine nijmegen log runs
    nijm_logs_missing = rbind(nijm_logs_missing, nijm_logs_newer_missing)
    nijm_logs_missing = nijm_logs_missing[!duplicated(nijm_logs_missing$Smrt),]
    dim(nijm_logs_missing)
    ggplot(nijm_logs_missing, aes(x=Date, y=Yield, color=Yield)) + geom_point(stat = 'identity', size=3) + ggtitle('Missing runs from Nijmegen only based on Log files') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
    # see how many runs with yield >600gb we have
    dim(nijm_logs_missing[which(nijm_logs_missing$Yield >= 600),])
    # combine all nijmegen runs
    all_nij_runs = rbind(nijm_logs, nijm_logs_newer)
    all_nij_runs = all_nij_runs[!duplicated(all_nij_runs$MovieID),]
    all_nij_runs$Site = 'Nijmegen'
    all_nij_runs$date = str_split_fixed(all_nij_runs$MovieID, '_', 3)[, 2]
    all_nij_runs$date = paste0('20', all_nij_runs$date)
    all_nij_runs$date = as.Date(all_nij_runs$date, "%Y%m%d")
    colnames(all_nij_runs) = c('Smrt', 'Yield', 'Site', 'Date')
    # check if we missed something
    tmp = fread('All_runs_Logs.txt', h=T)
    missing = tmp[which(!(tmp$MovieID %in% c(all_nij_runs$Smrt, all_run_logs$Smrt))),] # --> ALL is good here

    # maybe plot Yield over time (all smrt cells ever?)
    all_run_logs = fread('Run_logs_Ams.txt', h=T, stringsAsFactors=F, dec=',')
    colnames(all_run_logs) = c('MovieID', 'Yield')
    all_run_logs$Site = 'Amsterdam'
    all_run_logs$date = str_split_fixed(all_run_logs$MovieID, '_', 3)[, 2]
    all_run_logs$date = paste0('20', all_run_logs$date)
    all_run_logs$date = as.Date(all_run_logs$date, "%Y%m%d")
    colnames(all_run_logs) = c('Smrt', 'Yield', 'Site', 'Date')

    # combine amsterdam and nijmegen
    all_run_logs = rbind(all_run_logs, all_nij_runs)
    all_run_logs = all_run_logs[!duplicated(all_run_logs$Smrt),]
    all_run_logs$Yield = str_replace_all(all_run_logs$Yield, ',', '.')
    all_run_logs$Yield = as.numeric(all_run_logs$Yield)
    all_run_logs = all_run_logs[!is.na(all_run_logs$Yield),]
    all_run_logs = all_run_logs[order(all_run_logs$Date),]
    ggplot(all_run_logs, aes(x=Date, y=Yield, color=Site)) + geom_point(stat = 'identity', size=3, alpha = 0.5) + geom_smooth(method = 'loess') + ggtitle('SMRT Yield over time and sequencing centers') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
    
    #################################### -- Final integrating all things
    # AIM IS TO GENERATE AN OVERVIEW AT THE SAMPLE LEVEL
    # This is based on the merged samples
    merged_sequenced
    # This is based on the single smrt cells
    table_occurrences_filter
    # check whether we expected to sequence these
    tofix = table_occurrences_filter[which(!(table_occurrences_filter$Var1 %in% all_samples_expected$ID_GWAS)),]
    # 2 centenarians were (apparently) not supposed to be sequenced (nor sure why, nothing weird. maybe sample swap??), and 1 AD (EMC_NEURO study, probably also Anke/Annemieke). 3 Anke's samples are OK.
    table_occurrences_filter_info = merge(table_occurrences_filter, pheno_final_raw, by.x = 'Var1', by.y = 'ID_GWAS')
    table_occurrences_filter_info2 = merge(table_occurrences_filter_info, all_samples_expected, by.x = 'Var1', by.y = 'ID_GWAS', all.x=T)
    
    # then we need to go to Nijmegen data
    new_nijm = read.xlsx('Nijmegen_Sample_ID_and_RunIDs_Overview.xlsx')
    # exclude samples from the blood-brain-child project
    new_nijm = new_nijm[which(new_nijm$Project != 'BBC'),]
    # from this new list, exclude smrt cells already processed (in table occurrences)
    new_nijm$Done = NA
    for (i in 1:nrow(new_nijm)){
      # get run id
      tmp_smrt = new_nijm$MovieID[i]
      # check if was merged
      was_merged = merged_sequenced[grep(tmp_smrt, merged_sequenced$SMRT_CELLS),]
      was_proces = table_occurrences_filter[grep(tmp_smrt, table_occurrences_filter$Smrts),]
      if (nrow(was_merged) >0){
        new_nijm$Done[i] = 'merged'
      } else if (nrow(was_proces) >0){
        new_nijm$Done[i] = 'proces'
      } else{
        new_nijm$Done[i] = 'todo'
      }
    }
    table(new_nijm$Done, exclude=F)
    # then take only todo
    todo = new_nijm[which(new_nijm$Done == 'todo'),]
    todo$Sample = stringr::str_replace_all(todo$Sample, '-02', '')
    todo$Sample = stringr::str_replace_all(todo$Sample, '_19kb_30-07-2021', '')
    
    # read run stats
    run_stats = read.xlsx('Copy_nijmegen_runs_from_okt.xlsx', sheet = 2)
    table(todo$MovieID %in% run_stats$context)
    not_matching = todo[which(!todo$MovieID %in% run_stats$context),]
    todo = merge(todo, run_stats, by.x = 'MovieID', by.y = 'context', all.x=T)
    # 16 do not match between the run statistics from Nijmegen: this can be because the run statistics (run_stats) are only from october onwards, while the other file (todo) also includes other months, so that's ok
    table(todo$Sample %in% plates$barcode)
    # 7 samples sequenced were not in the plates -- maybe some change we did afterwards?
    mismatch = todo[which(!(todo$Sample %in% plates$barcode)),]
    # solve these 7
    missing1 = c('KG-011524', '17R1457', 'FR34137283', '202204011P1V', 'D01')
    missing2 = c('ND-001200', '11D01145', 'FR34137218', '202204011P1V', 'A01')
    missing3 = c('ND-000882', '10D04356', 'FR34137285', '202204011P1V', 'G04')
    missing4 = c('ND-000025', '09D04551', 'FR34137227', '202204011P1V', 'B06')
    missing5 = c('KG-011559', '17R1369', 'FR34137302', '202204011P1V', 'H01')
    df_missing = data.frame(plate = c(missing1[4], missing2[4], missing3[4], missing4[4], missing5[4]), biobank = c(missing1[2], missing2[2], missing3[2], missing4[2], missing5[2]), barcode = c(missing1[3], missing2[3], missing3[3], missing4[3], missing5[3]), kg_number = c(missing1[1], missing2[1], missing3[1], missing4[1], missing5[1]), plate_id = rep('other', 5))
    plates = rbind(plates, df_missing, use.names=TRUE)
    table(todo$Sample %in% plates$barcode)
    # merge with plates
    todo_pheno = merge(todo, plates, by.x = 'Sample', by.y = 'barcode', all.x=T)
    # add phenotypes
    table(todo_pheno$kg_number %in% all_samples_expected$"KG-number")
    todo_pheno_diagnosis = merge(todo_pheno, all_samples_expected, by.x = 'kg_number', by.y = 'KG-number')
    # merge samples with multiple smrt cells
    todo_pheno_diagnosis_merged = data.frame()
    for (i in 1:nrow(todo_pheno_diagnosis)){
      tmp = todo_pheno_diagnosis[which(todo_pheno_diagnosis$Sample == todo_pheno_diagnosis$Sample[i]),]
      tmp = tmp[!duplicated(tmp$MovieID),]
      tmp_df = data.frame(kg_number = unique(tmp$kg_number), sample = unique(tmp$Sample), Smrt = paste(tmp$MovieID, collapse=','), yield = sum(tmp$"Total.Bases.(Gb)"), id_gwas = unique(tmp$ID_GWAS), pheno = unique(tmp$pheno))
      todo_pheno_diagnosis_merged = rbind(todo_pheno_diagnosis_merged, tmp_df)
    }
    todo_pheno_diagnosis_merged = todo_pheno_diagnosis_merged[!duplicated(todo_pheno_diagnosis_merged$sample),]
    
    ######################################## 
    # SO NOW WE NEED TO PUT THINGS TOGETHER
    df_ready = merged_sequenced[, c('SMRT_CELLS', 'SAMPLE', 'COMBINED_COVERAGE', 'DIAGNOSIS', 'KG-number')]
    bbc = df_ready[is.na(df_ready$SMRT_CELLS),]; df_ready = df_ready[!is.na(df_ready$SMRT_CELLS),]
    colnames(df_ready) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
    df_single = table_occurrences_filter_info2[, c('Smrts', 'Var1', 'Cov', 'Pheno', 'KG-number')]          # modify here in case something is weird
    colnames(df_single) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
    df_nijmegen = todo_pheno_diagnosis[, c('MovieID', 'ID_GWAS', 'Total.Bases.(Gb)', 'pheno', 'kg_number')]
    colnames(df_nijmegen) = c('smrt', 'sample', 'yield', 'phenotype', 'kg_number')

    # would be nice to add the run yield for amsterdam runs
    ams_yield = data.table::fread('Run_logs_Ams.txt', h=T, stringsAsFactors=F, dec=',')
    colnames(ams_yield) = c('smrt', 'yield')
    dim(ams_yield)
    ams_yield = ams_yield[which(ams_yield$smrt != ''),]
    ams_yield$yield = as.numeric(stringr::str_replace_all(ams_yield$yield, ',', '.'))
    
    # nijmegen needs to be combined
    nijm_logs_cp = all_nij_runs
    colnames(nijm_logs_cp) = c('smrt', 'yield', 'site', 'date')
    run_stat_nj = run_stats[, c('context', 'Total.Bases.(Gb)')]
    run_stat_nj = run_stat_nj[which(!run_stat_nj$context %in% nijm_logs_cp$smrt),]
    colnames(run_stat_nj) = c('smrt', 'yield')
    nijm_logs_cp = rbind(nijm_logs_cp, run_stat_nj, fill=TRUE)

    # add yield for df_ready     
    # add a missing one
    nijm_logs_cp$site = NULL; nijm_logs_cp$date = NULL
    nijm_logs_cp = rbind(nijm_logs_cp, data.frame(smrt = c('m64102e_221115_215153', 'm64102e_221113_013123', 'm64102e_221111_143527', 'm64037e_221029_002710'), yield = c(640.103, 46.344, 85.575, 46.606)))
    df_ready$yield = NA
    for (i in 1:nrow(df_ready)){
      tmp = unlist(strsplit(df_ready$smrt[i], ','))
      # search in amsterdam
      tmp_yield = ams_yield[which(ams_yield$smrt %in% tmp),]
      # search in nijmegen
      tmp_yield2 = nijm_logs_cp[which(nijm_logs_cp$smrt %in% tmp),]
      tmp_yield2 = tmp_yield2[!duplicated(tmp_yield2$smrt),]
      # combine
      tmp_combined = rbind(tmp_yield, tmp_yield2)
      # check
      if (nrow(tmp_combined) != length(tmp)){
        print(paste('Check', i))
      } else {
        df_ready$yield[i] = sum(tmp_combined$yield)
      }
    }
    # add yield for processed
    df_single$yield = NA
    for (i in 1:nrow(df_single)){
      tmp = unlist(strsplit(df_single$smrt[i], ','))
      # search in amsterdam
      tmp_yield = ams_yield[which(ams_yield$smrt %in% tmp),]
      # search in nijmegen
      tmp_yield2 = nijm_logs_cp[which(nijm_logs_cp$smrt %in% tmp),]
      tmp_yield2 = tmp_yield2[!duplicated(tmp_yield2$smrt),]
      # combine
      tmp_combined = rbind(tmp_yield, tmp_yield2)
      # check
      if (nrow(tmp_combined) != length(tmp)){
        print(paste('Check', i))
      } else {
        df_single$yield[i] = sum(tmp_combined$yield)
      }
    }
    
    # combine all datasets
    df_ready$type = 'ready_merged'
    df_single$type = 'to_merge_or_single'
    df_nijmegen$type = 'nijmegen_to_run'
    df_combined_with_coverage = rbind(df_ready, df_single, use.names=TRUE)
    
    # show correlation between yield and coverage
    ggplot(data = df_combined_with_coverage, aes(x = coverage, y = yield)) + geom_point(stat = 'identity', size = 3) + geom_smooth() + xlab('Coverage of HiFi data') + ylab('Total Yield (GB)')

    # can we predict coverage of the nijmegen ones?
    model = lm(coverage ~ yield, data = df_combined_with_coverage)
    # fix some yield missing in df_nijmegen
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220922_051859')] = 647.275
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220624_092544')] = 579.398
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220621_142045')] = 550.855
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220816_194859.')] = 572.151
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_221002_025718')] = 34
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220625_185740')] = 458.482
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64102e_220622_234033')] = 583.096
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220627_144641')] = 527.324
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220829_124511')] = 408.805
    df_nijmegen$yield[which(df_nijmegen$smrt == 'm64037e_220830_222040')] = 440.76
    df_nijmegen$coverage = predict(model, newdata = df_nijmegen)
    
    # combined all together
    df_combined_all = rbind(df_combined_with_coverage, df_nijmegen)
    dim(df_combined_all)
    # we need to merge runs of the same samples again
    df_cleaned = data.frame()
    for (i in 1:nrow(df_combined_all)){
      if (!(df_combined_all$sample[i] %in% df_cleaned$sample)){
        tmp = df_combined_all[which(df_combined_all$sample == df_combined_all$sample[i]),]
        combined = data.frame(smrt = paste(tmp$smrt, collapse=','), sample = unique(tmp$sample), coverage = sum(tmp$coverage, na.rm=T), phenotype = unique(tmp$phenotype), kg_number = unique(tmp$kg_number), yield = sum(tmp$yield, na.rm=T), type = paste(tmp$type, collapse=','))
        df_cleaned = rbind(df_cleaned, combined)
      }
    }
    table(df_cleaned$phenotype)
    df_cleaned = df_cleaned[which(df_cleaned$smrt != 'TRUE'),]
    table(df_cleaned$phenotype)
    
    # plots
    # all samples, Amsterdam and Nijmegen together
    pdf('all_samples_with_prediction.pdf', height=7, width=10)
    ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
    dev.off()

    # same plot but color good and bad samples
    df_cleaned$Status = ifelse(df_cleaned$coverage >=12, 'Good', 'Re-do')
    pdf('all_samples_with_prediction_status.pdf', height=7, width=10)
    ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2, aes(col = Status)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
    dev.off()

    # low coverage
    # <12
    pdf('redo_with_predicted.pdf', height = 7, width = 10)
    low_cov = df_cleaned[which(df_cleaned$coverage <12),]
    table(low_cov$phenotype)
    ggplot(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
    dev.off()
    # also histogram
    pdf('redo_with_predicted_hist.pdf', height = 7, width = 10)
    ggplot(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = coverage, fill = phenotype)) + geom_histogram() + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Coverage of HiFi data') + ylab('Samples count')
    dev.off()
    save.image('20230512.RData')

    # output the list of samples with low coverage
    save.image('20230616.RData')
    LATEST_df_ready = df_cleaned

    # compare with previous workspace
    load('../data_20230512/20230512.RData')
    dim(df_cleaned)
    dim(df_cleaned[!duplicated(df_cleaned$sample),])
    dim(LATEST_df_ready)
    dim(LATEST_df_ready[!duplicated(LATEST_df_ready$sample),])
    diff = LATEST_df_ready[which(!(LATEST_df_ready$sample %in% df_cleaned$sample)),]
    head(diff)
###############################################################
###############################################################
# Updated on 30/06
    # MAKE AN OVERVIEW OF THE SEQUENCING STATUS
    # LIBRARIES
        library(data.table)
        library(plyr)
        library(openxlsx)
        library(stringr)
        library(ggplot2)

    # READ ALL SAMPLES WE WERE SUPPOSED TO DO
        all_chc = fread('chc_for_pacbio_withSCORE_and_KGnumb.txt', h=T)
        all_ad = fread('20200724_final_list_ADcases_completed.txt', h=T)
        # COMBINE THEM
        colnames(all_chc) = c('ID_100plus', 'ID_GWAS', 'pheno', 'gender', 'age', 'pca_outliers_gwas', 'pbmc', 'children', 'children_partners', 'siblings', 'trios', 'brain', 'score', 'KG-number')
        colnames(all_ad) = c('I_ID', 'pheno', 'info', 'gender', 'brain_year', 'apoe', 'age_death', 'age_diag', 'last_age_available', 'brain', 'KG-number')
        # for AD cases, need to add the GWAS ID
        load('/project/holstegelab/Share/gwas_array/mapping_files/phenotypes_20211027_All.Rdata')
        tmp_pheno = pheno_final_raw[, c('I_ID', 'ID_GWAS')]
        all_ad = merge(all_ad, tmp_pheno, by = 'I_ID', all.x = T)
        # there are 3 duplicates: that is, 3 samples (AD) with 2 ID_GWAS. 2 of them make sense, as they were genotyped multiple times, 1 is more controversial as it's from the same genotyping run but different phenotypes (male/female, I_ID = 2100). Keep all duplicates in for now.
        # merge
        all_samples_expected = rbind.fill(all_chc, all_ad)
        table(all_samples_expected$pheno)
        # we started of with 353 AD cases (3 duplicates) and 356 centenarians
    
    # READ ALL PLATES
        plates = data.frame()
        # plate 1 is in a different data, do it separately
        plate1 = read.xlsx('Plate1_QC.xlsx', sheet = 1)
        # restrict to plate, biobank, barcode and KG-number
        plate1 = plate1[, c('PlateIDtarget', 'SampleID', 'Barcode', 'Individual.ID')]
        # rename columns
        colnames(plate1) = c('plate', 'biobank', 'barcode', 'kg_number')
        # add plate id
        plate1$plate_id = 1
        # add to plates
        plates = rbind(plates, plate1)
        for (p in 2:7){
            # read sheet
            tmp = read.xlsx('All_plates_QC.xlsx', sheet = p)
            # restrict to plate, biobank, barcode and KG-number
            tmp = tmp[, c('Plate', 'Biobank.fractie', 'BC', 'KG-number')]
            # rename columns
            colnames(tmp) = c('plate', 'biobank', 'barcode', 'kg_number')
            # add plate id
            tmp$plate_id = p
            # add to plates
            plates = rbind(plates, tmp)
        }
        table(plates$plate_id)
        # the plates contain 653 samples in total
    
    # NOW READ ALL SEQUENCED SAMPLES BASED ON THE COVERAGE SUMMARY FILE
        # combine all coverage information
        coverage_flist = system('ls /project/holstegelab/Share/pacbio/data_processed/coverage_freezes/', intern=T)
        coverage_combined = data.frame()
        for (f in coverage_flist){
            tmp = fread(paste0('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/', f), dec=',')
            coverage_combined = rbind(coverage_combined, tmp, fill=T)
        }
        # check and exclude duplicates
        coverage_combined = coverage_combined[order(coverage_combined$GLOBAL_MEDIAN_READ_LENGTH),]
        coverage_combined = coverage_combined[!duplicated(coverage_combined$movie_id),]
        coverage_combined = coverage_combined[which(coverage_combined$movie_id != TRUE),]
        # read the merged samples
        merged_sequenced = fread('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', h=T, stringsAsFactors=F, sep="\t", dec = ',')
        merged_sequenced = merged_sequenced[!is.na(merged_sequenced$COMBINED_COVERAGE),]
        merged_sequenced = merged_sequenced[order(-merged_sequenced$COMBINED_COVERAGE),]
        merged_sequenced = merged_sequenced[!duplicated(merged_sequenced$SAMPLE),]
        # also read blood-brain-child samples
        bbc_samples = fread('/project/holstegelab/Software/snakemake_pipeline/sample_check_data/20210615_overview_genetics.txt', h=T, stringsAsFactors=F)
        bbc_samples = bbc_samples[!is.na(bbc_samples$ID_PACBIO),]
        bbc_samples = bbc_samples[, c('ID_GWAS')]
        bbc_samples$DIAGNOSIS = 'Blood_Brain_Child'
        colnames(bbc_samples) = c('SAMPLE', 'DIAGNOSIS')
        merged_sequenced = rbind.fill(merged_sequenced, bbc_samples)
        table(merged_sequenced$DIAGNOSIS)
        # in total, ready for analysis there are 91 AD cases + 114 Centenarians + 10 more centenarians (blood-brain-child project) + 10 children (blood-brain-child project)
    
    # NOW LABEL THE SAMPLES SEQUENCED IN ALL_SAMPLES WE WERE SUPPOSED TO SEQUENCE
        table(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)
        # 214/225 are OK -- let's check the remaining 11
        non_matching = merged_sequenced[which(!(merged_sequenced$SAMPLE %in% all_samples_expected$ID_GWAS)),]
        # add phenotypes
        non_matching_pheno = merge(non_matching, pheno_final_raw, by.x = 'SAMPLE', by.y = 'ID_GWAS')
        # it's ok -- these 10 are the children of the centenarians
        # one is a centenarian that apparently we were not supposed to sequence but we did -- peace and love
        # So all samples we sequenced until now were in the list of samples to be sequenced
    
    # NOW CHECK WITH RESPECT TO THE PLATES
        # add kg-number to the sequenced samples
        kg_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
        merged_sequenced = merge(merged_sequenced, kg_info, by.x = 'SAMPLE', by.y = 'ID_GWAS', all.x = T)
        # fix the missing KG
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '8988312')] = 'KG-012463'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '17R2345')] = 'KG-012023'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '17R1919')] = 'KG-011789'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '17R2183')] = 'KG-011940'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '18R1671')] = 'KG-013139'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '18R1705')] = 'KG-013154'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '18R1739')] = 'KG-013176'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '18R2415')] = 'KG-013533'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '19R0074')] = 'KG-013617'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '19R0807')] = 'KG-013994'
        merged_sequenced$"KG-number"[which(merged_sequenced$SAMPLE == '8629102')] = 'KG-012617'
        table(merged_sequenced$DIAGNOSIS)
        # merge sequenced samples with plate information
        merged_sequenced_plate = merge(merged_sequenced, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
        # add phenotypes
        tmp_pheno = pheno_final_raw[, c('ID_GWAS', 'diagnosis', 'ID_100plus', 'I_ID')]
        merged_sequenced_plate = merge(merged_sequenced_plate, tmp_pheno, by.x = 'SAMPLE', by.y = 'ID_GWAS')
        table(merged_sequenced_plate$DIAGNOSIS)
        # there are duplicates, we should clean them up
        dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
        merged_sequenced_plate = merged_sequenced_plate[order(-merged_sequenced_plate$COMBINED_COVERAGE),]
        merged_sequenced_plate = merged_sequenced_plate[!duplicated(merged_sequenced_plate$SAMPLE),] 
        dups = merged_sequenced_plate[duplicated(merged_sequenced_plate$SAMPLE),]
        # find samples with missing plate id
        merged_sequenced_plate[is.na(merged_sequenced_plate$plate_id),]
        table(merged_sequenced_plate$plate_id, exclude=F)
        merged_sequenced_plate$SMRT_CELL_N = factor(merged_sequenced_plate$SMRT_CELL_N)
        # 32 in total: 20 are centenarians from the blood-brain-child project (centenarians + childre), 4 are additional centenarians + 8 AD
    
    # PLOT THE COMBINED_COVERAGE OF SAMPLES
        pdf('coverage_merged_hq_samples.pdf', height = 7, width = 10)
        ggplot(merged_sequenced_plate[which(merged_sequenced_plate$DIAGNOSIS != 'Blood_Brain_Child'),], aes(x = DIAGNOSIS, y = COMBINED_COVERAGE, fill = DIAGNOSIS)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(shape = SMRT_CELL_N), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Combined coverage of HiFi reads')
        dev.off()
        # coverage stats
        mean(merged_sequenced_plate$COMBINED_COVERAGE[which(merged_sequenced_plate$DIAGNOSIS == 'Probable_AD')])
        sd(merged_sequenced_plate$COMBINED_COVERAGE[which(merged_sequenced_plate$DIAGNOSIS == 'Probable_AD')])
        mean(merged_sequenced_plate$COMBINED_COVERAGE[which(merged_sequenced_plate$DIAGNOSIS == 'Centenarian')])
        sd(merged_sequenced_plate$COMBINED_COVERAGE[which(merged_sequenced_plate$DIAGNOSIS == 'Centenarian')])
    
    # LOOK AT SMRT CELLS >600GB YIELD
        dim(coverage_combined)
        # 757 good smrt cells sequenced thus far --> make output to update my Excel
        write.table(coverage_combined, '/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/2023-06-30_combined_coverage_all.txt', quote=F, row.names=F, sep='\t', dec=',')
        # ok now exclude samples combined (from merged_sequenced_plate)
        smrt_coverage_all_notcombined = coverage_combined[which(!(coverage_combined$ID_GWAS %in% merged_sequenced_plate$SAMPLE)),]
        dim(smrt_coverage_all_notcombined)
        smrt_coverage_all_notcombined = smrt_coverage_all_notcombined[which(!(smrt_coverage_all_notcombined$ID_GWAS %in% bbc_samples$SAMPLE)),]
        table_occurrences = data.frame(table(smrt_coverage_all_notcombined$ID_GWAS))
        table_occurrences = table_occurrences[order(-table_occurrences$Freq),]
        head(table_occurrences)
        # calculate total coverage of these
        table_occurrences$Cov = NA
        table_occurrences$Pheno = NA
        table_occurrences$Smrts = NA
        for (i in 1:nrow(table_occurrences)){
            tmp = smrt_coverage_all_notcombined[which(smrt_coverage_all_notcombined$ID_GWAS == table_occurrences$Var1[i]),]
            table_occurrences$Cov[i] = sum(tmp$GLOBAL_COVERAGE)
            table_occurrences$Pheno[i] = unique(tmp$diagnosis)
            table_occurrences$Smrts[i] = paste(tmp$movie_id, collapse = ',')
        }
        table_occurrences = table_occurrences[order(-table_occurrences$Cov),]
        head(table_occurrences, 10)

    # ALSO CHECK IN THE MOST UPDATED LIST OF COMBINED SAMPLES (THIS CAN ADD SOMETHING IN CASE SAMPLES ARE BEING MERGED AT THE MOMENT)
        samples_merged_chc = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/centenarian/*merged.hifi.chm13.bam', intern=T)
        samples_merged_ad = system('ls /project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/ad/*merged.hifi.chm13.bam', intern=T)
        df_chc = data.frame(sample = samples_merged_chc, pheno = 'chc')
        df_ad = data.frame(sample = samples_merged_ad, pheno = 'ad')
        df_combined = rbind(df_chc, df_ad)
        for (i in 1:nrow(df_combined)){
            tmp = str_replace_all(strsplit(df_combined$sample[i], '/')[[1]][length(strsplit(df_combined$sample[i], '/')[[1]])], '.merged.hifi.chm13.bam', '')
            df_combined$sample[i] = tmp
        }
        table(df_combined$pheno)
        # exclude from the table of non combined any of these actually combined
        dim(table_occurrences)
        table_occurrences_filter = table_occurrences[which(!(table_occurrences$Var1 %in% df_combined$sample)), ]
        dim(table_occurrences_filter)

    # TAKE OUT THE SAMPLES WITH COVERAGE >12
        coverage_enouth = table_occurrences_filter[which(table_occurrences_filter$Cov >= 12), ]
        dim(coverage_enouth)
        coverage_enouth$Freq = factor(coverage_enouth$Freq)
        coverage_enouth$Pheno[which(coverage_enouth$Pheno == 'Control_100plus')] = 'Centenarian'
        coverage_enouth$Pheno[which(coverage_enouth$Pheno == 'MCI')] = 'Probable_AD'
    
    # PLOT NON-MERGED SAMPLES (BOTH ENOUGH COVERAGE AND TO RESEQUENCE)
        pdf('all_samples_with_data.pdf', height = 7, width = 10)
        table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno == 'Control_100plus')] = 'Centenarian'
        table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno %in% c('AD_othersite', 'MCI'))] = 'Probable_AD'
        table_occurrences_filter$Pheno[which(table_occurrences_filter$Pheno %in% c('NIID', 'Dementie_anders'))] = 'Anke'
        table_occurrences_filter = table_occurrences_filter[which(table_occurrences_filter$Pheno != 'family_100plus'),]
        table_occurrences_filter$Freq = factor(table_occurrences_filter$Freq)
        ggplot(table_occurrences_filter[which(table_occurrences_filter$Pheno %in% c('Centenarian', 'Probable_AD')),], aes(x = Pheno, y = Cov, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', aes(color = Freq), width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + ylab('Coverage of HiFi reads') + xlab('Phenotype') + ylim(0, 45) + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
        dev.off()
    # OVERLAP WITH PLATES
        plates_info = all_samples_expected[, c('ID_GWAS', 'KG-number')]
        coverage_enouth = merge(coverage_enouth, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
        coverage_enouth_plate = merge(coverage_enouth, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
        # check plates
        table(merged_sequenced_plate$plate_id, exclude=F)
        table(coverage_enouth_plate$plate_id, exclude=F)
        # all samples with data
        table_occurrences_filter_info = merge(table_occurrences_filter, plates_info, by.x = 'Var1', by.y = 'ID_GWAS', all.x = T)
        table_occurrences_filter_info_plate = merge(table_occurrences_filter_info, plates, by.x = 'KG-number', by.y = 'kg_number', all.x = T)
        table(table_occurrences_filter_info_plate$plate_id, exclude=F)

    # NEXT, GO TO NIJMEGEN STATS
        # read older logs
        nijm_logs = fread('Nijmegen_Log.txt', h=T, stringsAsFactors=F, dec=',')
        colnames(nijm_logs) = c('MovieID', 'Yield')
        # add the very new data
        nijm_logs_newer = fread('20230630_Nijmegen_Logs.txt', h=T, stringsAsFactors=F, dec=',')
        dim(nijm_logs_newer)
        # check if there's something in old not in new
        table(nijm_logs$MovieID %in% nijm_logs_newer$MovieID)
        # only 1 missing --> combine and remove duplicates
        nijm_logs_all = rbind(nijm_logs_newer, nijm_logs, fill=T)
        nijm_logs_all = nijm_logs_all[!duplicated(nijm_logs_all$MovieID),]
        dim(nijm_logs_all)
        # exclude runs we have processed
        nijm_logs_all_missing = nijm_logs_all[which(!(nijm_logs_all$MovieID %in% coverage_combined$movie_id)),]
        # for some reason some runs are not in the coverage_combined data
        merged_samples_smrtID = unlist(strsplit(merged_sequenced$SMRT_CELLS, ','))
        nijm_logs_all_missing = nijm_logs_all_missing[which(!(nijm_logs_all_missing$MovieID %in% merged_samples_smrtID)),]
        dim(nijm_logs_all_missing)
        # n=95 (1 run from june, few from October/November and from January until now)
        # plot yield and date
        nijm_logs_all_missing$date = str_split_fixed(nijm_logs_all_missing$MovieID, '_', 3)[, 2]
        nijm_logs_all_missing$date = paste0('20', nijm_logs_all_missing$date)
        nijm_logs_all_missing$date = as.Date(nijm_logs_all_missing$date, "%Y%m%d")
        colnames(nijm_logs_all_missing) = c('Smrt', 'Movie', 'Yield', 'Sample', 'Date')
        # check which samples these are
        nijm_logs_all_missing_plates = merge(nijm_logs_all_missing, plates, by.x = 'Sample', by.y = 'barcode', all.x = T)
        # check mismatches
        nijm_logs_all_missing_plates[is.na(nijm_logs_all_missing_plates$biobank),]
        nijm_logs_all_missing_plates$kg_number[which(nijm_logs_all_missing_plates$Sample == 'FR34137226')] = 'KG-008407'
        nijm_logs_all_missing_plates$biobank[which(nijm_logs_all_missing_plates$Sample == 'FR34137226')] = '15R3339'
        #
        nijm_logs_all_missing_plates$kg_number[which(nijm_logs_all_missing_plates$Sample == 'FR34137227')] = 'ND-000025'
        nijm_logs_all_missing_plates$biobank[which(nijm_logs_all_missing_plates$Sample == 'FR34137227')] = '09D04551'
        #
        nijm_logs_all_missing_plates$kg_number[which(nijm_logs_all_missing_plates$Sample == 'FR34137228')] = 'ND-000877'
        nijm_logs_all_missing_plates$biobank[which(nijm_logs_all_missing_plates$Sample == 'FR34137228')] = '10D04351'
        #
        nijm_logs_all_missing_plates$kg_number[which(nijm_logs_all_missing_plates$Sample == 'FR34137283')] = 'KG-011524'
        nijm_logs_all_missing_plates$biobank[which(nijm_logs_all_missing_plates$Sample == 'FR34137283')] = '17R1457'
        #
        nijm_logs_all_missing_plates$kg_number[which(nijm_logs_all_missing_plates$Sample == 'FR34137291')] = 'ND-000044'
        nijm_logs_all_missing_plates$biobank[which(nijm_logs_all_missing_plates$Sample == 'FR34137291')] = '09D04570'
        #
        nijm_logs_all_missing_plates$kg_number[which(nijm_logs_all_missing_plates$Sample == 'FR34137302')] = 'KG-011559'
        nijm_logs_all_missing_plates$biobank[which(nijm_logs_all_missing_plates$Sample == 'FR34137302')] = '17R1369'
        nijm_logs_all_missing_plates[is.na(nijm_logs_all_missing_plates$biobank),]
        # add phenotype
        table(nijm_logs_all_missing_plates$kg_number %in% c(all_samples_expected$"KG-number"))
        nijm_logs_all_missing_plates = merge(nijm_logs_all_missing_plates, all_samples_expected, by.x = 'kg_number', by.y = 'KG-number')
        # plot
        ggplot(nijm_logs_all_missing_plates, aes(x=Date, y=Yield, color=Yield)) + geom_point(stat = 'identity', size=3, aes(shape = pheno)) + ggtitle('Missing runs from Nijmegen only based on Log files') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
        # see how many runs with yield >600gb we have
        dim(nijm_logs_all_missing_plates[which(nijm_logs_all_missing_plates$Yield >= 600),])

    # PLOT YIELD OVER TIME
        all_run_logs = fread('20230630_Ams_Logs.txt', h=T, stringsAsFactors=F, dec=',')
        colnames(all_run_logs) = c('MovieID', 'Yield')
        all_run_logs$Site = 'Amsterdam'
        nijm_logs_all$Site = 'Nijmegen'
        # combine amsterdam and nijmegen
        all_run_logs = rbind(all_run_logs, nijm_logs_all, fill=T)
        all_run_logs$date = str_split_fixed(all_run_logs$MovieID, '_', 3)[, 2]
        all_run_logs$date = paste0('20', all_run_logs$date)
        all_run_logs$date = as.Date(all_run_logs$date, "%Y%m%d")
        all_run_logs = all_run_logs[!duplicated(all_run_logs$MovieID),]
        all_run_logs$Yield = str_replace_all(all_run_logs$Yield, ',', '.')
        all_run_logs$Yield = as.numeric(all_run_logs$Yield)
        all_run_logs = all_run_logs[!is.na(all_run_logs$Yield),]
        all_run_logs = all_run_logs[order(all_run_logs$date),]
        ggplot(all_run_logs, aes(x=date, y=Yield, color=Site)) + geom_point(stat = 'identity', size=3, alpha = 0.5) + geom_smooth(method = 'loess') + ggtitle('SMRT Yield over time and sequencing centers') + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16))
    
    # DATA PRODUCED THUS FAR IN THIS OVERVIEW
        # merged samples
        merged_sequenced
        # non-merged samples
        table_occurrences_filter
        # missing from nijmegen
        nijm_logs_all_missing_plates
        # all from nijmegen
        nijm_logs_all
        # all logs
        all_run_logs
    
    # PUT THINGS TOGETHER
        df_ready = merged_sequenced[, c('SMRT_CELLS', 'SAMPLE', 'COMBINED_COVERAGE', 'DIAGNOSIS', 'KG-number')]         # samples ready for analyses
        bbc = df_ready[is.na(df_ready$SMRT_CELLS),]; df_ready = df_ready[!is.na(df_ready$SMRT_CELLS),]
        colnames(df_ready) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
        df_single = table_occurrences_filter_info[, c('Smrts', 'Var1', 'Cov', 'Pheno', 'KG-number')]          # modify here in case something is weird
        colnames(df_single) = c('smrt', 'sample', 'coverage', 'phenotype', 'kg_number')
        df_nijmegen = nijm_logs_all_missing_plates[, c('Movie', 'ID_GWAS', 'Yield', 'pheno', 'kg_number')]
        colnames(df_nijmegen) = c('smrt', 'sample', 'yield', 'phenotype', 'kg_number')
        # extract amsterdam yields
        ams_yield = all_run_logs[which(all_run_logs$Site == 'Amsterdam')]
        # merge samples with multiple smrt cells from nijmegen
        nijm_logs_all_missing_plates_sample = data.frame()
        for (s in unique(nijm_logs_all_missing_plates$ID_GWAS)){
            tmp = nijm_logs_all_missing_plates[which(nijm_logs_all_missing_plates$ID_GWAS == s),]
            tmp = tmp[!duplicated(tmp$Movie),]
            tmp_df = data.frame(kg_number = unique(tmp$kg_number), sample = unique(tmp$Sample), biobank = unique(tmp$biobank), Smrt = paste(tmp$Movie, collapse=','), yield = sum(tmp$Yield), id_gwas = unique(tmp$ID_GWAS), pheno = unique(tmp$pheno))
            nijm_logs_all_missing_plates_sample = rbind(nijm_logs_all_missing_plates_sample, tmp_df)
        }
        nijm_logs_all_missing_plates_sample[duplicated(nijm_logs_all_missing_plates_sample$sample),]
        df_nijmegen = nijm_logs_all_missing_plates_sample[, c('Smrt', 'id_gwas', 'yield', 'pheno', 'kg_number')]; colnames(df_nijmegen) = c('smrt', 'sample', 'yield', 'phenotype', 'kg_number')

        # add the run yields to the amsterdam runs to merge and the merged
        df_ready$yield = NA
        for (i in 1:nrow(df_ready)){
            tmp = unlist(strsplit(df_ready$smrt[i], ','))
            # search yields
            tmp_yield = all_run_logs[which(all_run_logs$MovieID %in% tmp),]
            # check
            if (nrow(tmp_yield) != length(tmp)){
                print(paste('Check', i))
            } else {
                df_ready$yield[i] = sum(tmp_yield$Yield)
            }
        }
        df_ready_qc = df_ready[!is.na(df_ready$yield),]
        # add yield for processed
        df_single$yield = NA
        for (i in 1:nrow(df_single)){
            tmp = unlist(strsplit(df_single$smrt[i], ','))
            tmp_yield = all_run_logs[which(all_run_logs$MovieID %in% tmp),]
            # check
            if (nrow(tmp_yield) != length(tmp)){
                print(paste('Check', i))
            } else {
                df_single$yield[i] = sum(tmp_yield$Yield)
            }
        }
    
        # combine all datasets
        df_ready$type = 'ready_merged'
        df_single$type = 'to_merge_or_single'
        df_nijmegen$type = 'nijmegen_to_run'
        df_combined_with_coverage = rbind(df_ready, df_single)
    
        # show correlation between yield and coverage
        ggplot(data = df_combined_with_coverage, aes(x = coverage, y = yield)) + geom_point(stat = 'identity', size = 3) + geom_smooth() + xlab('Coverage of HiFi data') + ylab('Total Yield (GB)')

        # can we predict coverage of the nijmegen ones?
        model = lm(coverage ~ yield, data = df_combined_with_coverage)
        # fix some yield missing in df_nijmegen
        df_nijmegen$coverage = predict(model, newdata = df_nijmegen)
    
        # combined all together: merged, ready and nijmegen to arrive
        df_combined_all = rbind(df_combined_with_coverage, df_nijmegen)
        dim(df_combined_all)
        # we need to merge runs of the same samples again
        df_cleaned = data.frame()
        for (i in 1:nrow(df_combined_all)){
            if (!(df_combined_all$sample[i] %in% df_cleaned$sample)){
                tmp = df_combined_all[which(df_combined_all$sample == df_combined_all$sample[i]),]
                combined = data.frame(smrt = paste(tmp$smrt, collapse=','), sample = unique(tmp$sample), coverage = sum(tmp$coverage, na.rm=T), phenotype = unique(tmp$phenotype), kg_number = unique(tmp$kg_number), yield = sum(tmp$yield, na.rm=T), type = paste(tmp$type, collapse=','))
                df_cleaned = rbind(df_cleaned, combined)
            }
        }
        dim(df_cleaned)
        table(df_cleaned$phenotype)
        df_cleaned = df_cleaned[which(df_cleaned$smrt != 'TRUE'),]
        table(df_cleaned$phenotype)
    
    # LAST THING IS TO CHECK: WHICH DATA WERE ALREADY MERGED BUT WE NEED TO RE-DO BECAUSE OF A NEW SAMPLE
        # need to read all the config files
        all_config_files = system('ls /project/holstegelab/Software/snakemake_pipeline/config/config_merge/config*.yml', intern=T)
        # read all config files and store the results
        df = data.frame()
        for (f in all_config_files){
            tmp_f = read.table(f)
            smrt_cells = unlist(strsplit(tmp_f$V3[1], ','))
            smrt_cells_combined = ''
            n_smrt = 0
            # split smrt cells
            for (x in smrt_cells){
                n_smrt = n_smrt + 1
                smrt_cells_combined = ifelse(smrt_cells_combined == '', unlist(strsplit(smrt_cells[1], '/'))[length(unlist(strsplit(smrt_cells[1], '/')))], paste(smrt_cells_combined, unlist(strsplit(smrt_cells[1], '/'))[length(unlist(strsplit(smrt_cells[1], '/')))], sep=','))
            }
            # take sample name
            sample_name = unlist(strsplit(tmp_f$V3[2], '/'))[length(unlist(strsplit(tmp_f$V3[2], '/')))]
            # make dataframe
            df = rbind(df, data.frame(sample = sample_name, smrt_cells = smrt_cells_combined, n_smrt_cells = n_smrt))
        }

    # PLOT ALL SAMPLES FROM AMSTERDAM AND NIJMEGEN
        pdf('all_samples_with_prediction.pdf', height=7, width=10)
        ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
        dev.off()
        mean(df_cleaned$coverage[which(df_cleaned$phenotype == 'Centenarian')]); sd(df_cleaned$coverage[which(df_cleaned$phenotype == 'Centenarian')]); median(df_cleaned$coverage[which(df_cleaned$phenotype == 'Centenarian')])
        mean(df_cleaned$coverage[which(df_cleaned$phenotype == 'Probable_AD')]); sd(df_cleaned$coverage[which(df_cleaned$phenotype == 'Probable_AD')]); median(df_cleaned$coverage[which(df_cleaned$phenotype == 'Probable_AD')])

    # SAME PLOT BUT COLORING GOOD AND BAD SAMPLES
        df_cleaned$Status = ifelse(df_cleaned$coverage >=12, 'Good', 'Re-do')
        pdf('all_samples_with_prediction_status.pdf', height=7, width=10)
        ggplot(df_cleaned[which(df_cleaned$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2, aes(col = Status)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
        dev.off()
        table(df_cleaned$Status, df_cleaned$phenotype)

    # SAMPLES WITH LOW COVERAGE
        # <12
        pdf('redo_with_predicted.pdf', height = 7, width = 10)
        low_cov = df_cleaned[which(df_cleaned$coverage <12),]
        table(low_cov$phenotype)
        ggplot(low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),], aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 2) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
        dev.off()
        # color based on barcoding or entire resequencing
        low_cov_ad_chc = low_cov[which(low_cov$phenotype %in% c('Centenarian', 'Probable_AD')),]
        low_cov_ad_chc$Barcoding = 'No BC (Resequence)'
        low_cov_ad_chc$Barcoding[which(low_cov_ad_chc$coverage >=9)] = 'BC (2 samples)'
        low_cov_ad_chc$Barcoding[which(low_cov_ad_chc$coverage >=10)] = 'BC (3 samples)'
        low_cov_ad_chc$Barcoding[which(low_cov_ad_chc$coverage >=11.5)] = 'BC (4 samples)'
        # plot again
        ggplot(low_cov_ad_chc, aes(x = phenotype, y = coverage, fill = phenotype)) + geom_violin() + geom_boxplot(width=0.05, outlier.shape=NA) + geom_jitter(stat = 'identity', width = 0.25, size = 3, aes(color = Barcoding)) + theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"), legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=16)) + scale_fill_brewer(name = "Title", direction = -1, palette = "Set1") + xlab('Phenotype') + ylab('Coverage of HiFi data')
        table(low_cov_ad_chc$Barcoding)

    # OUTPUT LIST OF SAMPLES WITH LOW COVERAGE
        # check how many samples to re-do we have already done
        already_done = low_cov_ad_chc[which(low_cov_ad_chc$sample %in% df$sample),]
        low_cov_ad_chc = low_cov_ad_chc[which(!(low_cov_ad_chc$sample %in% already_done$sample)),]
    # FINALLY, SOME MANUAL SAMPLES TO CHECK (RUN IN AMSTERDAM, NO RESULTS YET, ONLY YIELD)
        # FR34156836
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156836')]
        df_cleaned$yield[which(df_cleaned$kg_number == tmp_kg)] # new run has 508 yield --> sample is OK, remove from list to do
        low_cov_ad_chc[which(low_cov_ad_chc$kg_number == tmp_kg),]
        low_cov_ad_chc$kg_number[is.na(low_cov_ad_chc$kg_number)] = 'NA'
        low_cov_ad_chc = low_cov_ad_chc[which(!(low_cov_ad_chc$kg_number == tmp_kg)),]
        # FR34156874
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156874')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # this was already merged with 17x coverage --> more data is guess
        low_cov_ad_chc[which(low_cov_ad_chc$kg_number == tmp_kg),]
        # FR34156826
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156826')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # new run has 386 yield --> this should be repeated entirely
        # FR34156833
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156833')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # new run has 471 yield --> this should be OK, remove from the list to do
        low_cov_ad_chc = low_cov_ad_chc[which(!(low_cov_ad_chc$kg_number == tmp_kg)),]
        dim(low_cov_ad_chc)
        # FR34111758
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34111758')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # new run has 619 yield --> this was not merged yet, but it is OK.
        # FR34156812
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156812')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # new run has 539 yield --> this should be OK now, remove from the list to do
        low_cov_ad_chc = low_cov_ad_chc[which(!(low_cov_ad_chc$kg_number == tmp_kg)),]
        dim(low_cov_ad_chc)
        # FR34156844
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156844')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # new run has 520 yield --> this should be OK now, remove from the list to do
        low_cov_ad_chc = low_cov_ad_chc[which(!(low_cov_ad_chc$kg_number == tmp_kg)),]
        dim(low_cov_ad_chc)
        # FR34156843
        tmp_kg = plates$kg_number[which(plates$barcode == 'FR34156843')]
        df_cleaned[which(df_cleaned$kg_number == tmp_kg),] # new run has 516 yield --> this should be OK now, remove from the list to do
        low_cov_ad_chc = low_cov_ad_chc[which(!(low_cov_ad_chc$kg_number == tmp_kg)),]
        dim(low_cov_ad_chc)
        # ok good, have a look also at the ones planned
        planned = c('FR34156870', 'FR34156866', 'FR34111683', 'FR34111776', 'FR34156802', 'FR34111714', 'FR34162548', 'FR34162530', 'FR34162546', 'FR34162557', 'FR34137245', 'FR34162512', 'FR34137301', 'FR34137250', 'FR34111764', 'FR34137229')
        tmp_kgs = plates[which(plates$barcode %in% planned),]
        df_cleaned[which(df_cleaned$kg_number %in% tmp_kgs$kg_number),]
        planned[which(!(planned %in% tmp_kgs$barcode))]
        # check missing
        # FR34137245 --> ND-001337
        df_cleaned[which(df_cleaned$kg_number == 'ND-001337'),]
        # FR34137301 --> ND-000230
        df_cleaned[which(df_cleaned$kg_number == 'ND-000230'),]
        # FR34137250 --> ND-001089
        df_cleaned[which(df_cleaned$kg_number == 'ND-001089'),]
        # FR34137229 --> ND-000996
        df_cleaned[which(df_cleaned$kg_number == 'ND-000996'),]
        # everything should be alright
        write.table(low_cov_ad_chc, '20230704_samples_to_resequence.txt', quote=F, row.names=F, sep="\t", dec=',')
    # SAVE IMAGE
        save.image('20230630.RData')
        


    # COMPARE WITH PREVIOUS LIST OF SAMPLES
        prev_list = read.table('previous_sample_set.txt', h=F)
        table(prev_list$V1 %in% low_cov_ad_chc$kg_number)
        diff = prev_list[which(!(prev_list$V1 %in% low_cov_ad_chc$kg_number)),]
        table(prev_list$V1 %in% df_cleaned$kg_number)
        df_cleaned[which(df_cleaned$kg_number %in% diff),]