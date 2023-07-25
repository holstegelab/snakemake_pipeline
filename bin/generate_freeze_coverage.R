# Script to generate freeze of coverage information
#
# Briefly, takes coverage stats, and merge with additional information ragarding sample ID, sample type, project and sequencing center

# Libraries
    library(data.table)
    library(stringr)

# Main
# 1. read input coverage stats and remove duplicates
    cat('## Reading coverage statistics..\n')
    coverage_stats = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_smrt_cells.txt', h=T, stringsAsFactors=F, sep='\t')
    coverage_stats = coverage_stats[!duplicated(coverage_stats),]

# 2. now that we are moving processed data to dcache, we need to adopt a different strategy -- probably good to store all sample-check files somewhere
    # list all sample-check-files in dcache
    dcache_config = '/project/holstegelab/Data/dcache_processed.conf'
    all_files_processed_dcache1 = system(paste0("rclone ls --config ", dcache_config, " dcache_processed:ccs/ad_centenarians/"), intern = T)
    all_files_processed_dcache2 = system(paste0("rclone ls --config ", dcache_config, " dcache_processed:ccs/blood_brain_child/"), intern = T)
    all_files_processed_dcache3 = system(paste0("rclone ls --config ", dcache_config, " dcache_processed:ccs/other/"), intern = T)
    all_files_processed_dcache = c(all_files_processed_dcache1, all_files_processed_dcache2, all_files_processed_dcache3)
    for (i in 1:length(all_files_processed_dcache)){
        tmp = strsplit(all_files_processed_dcache[i], ' ')[[1]][length(strsplit(all_files_processed_dcache[i], ' ')[[1]])]
        if (length(grep('sample', tmp)) >0){
            system(paste0('rclone copy -vv --progress --multi-thread-streams 1 --config ', dcache_config, ' dcache_processed:ccs/ad_centenarians/', tmp, ' /project/holstegelab/Share/pacbio/data_processed/sample_checks/'))
        }
    }
    # then there are also some data on disk
    files_disk_p1 = system('ls /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/*sample*', intern = T)
    files_disk_p2 = system('ls /project/holstegelab/Share/pacbio/data_processed/nijmegen/*sample*', intern = T)
    files_disk_p3 = system('ls /project/holstegelab/Share/pacbio/data_processed/other_samples/*sample*', intern = T)
    files_disk_p4 = system('ls /project/holstegelab/Share/pacbio/data_processed/Anke_samples/*sample*', intern = T)
    files_disk_p5 = system('find /project/holstegelab/Share/pacbio/data_processed/blood_brain_child/ -name "*sample*"', intern = T)
    all_files_disk = c(files_disk_p1, files_disk_p2, files_disk_p3, files_disk_p4, files_disk_p5)
    for (f in all_files_disk){ system(paste0('cp ', f, ' /project/holstegelab/Share/pacbio/data_processed/sample_checks/')) }

# 3. Open all sample match files that we will merge with the coverage data
    data_path = '/project/holstegelab/Share/pacbio/data_processed/sample_checks/'
    all_sample_checks = system(paste0('ls ', data_path), intern=T)
    all_checks_df = data.frame()
    for (f in all_sample_checks){
        # open file
        tmp = fread(paste0(data_path, f), h=T, stringsAsFactors=F)
        # sort by percentage of identity
        tmp = tmp[order(tmp$PERC_HOMOLOGY),]
        # take percentage higher than 95%
        likely_match = tmp[which(as.numeric(tmp$PERC_HOMOLOGY) >= 0.94),]
        # extract smrt id of the run
        smrt_id = str_replace_all(f, '.ccs.primrose.hifi.sample.txt', '')
        # check if there are no hits, and if there are more than 1 hit
        if (nrow(likely_match) == 0){
            print(paste0('!!! no hits for ', f, ': likely coverage is too low or the sample is unknown. Now trying to figure it out.'))
            # check what samples this was
            if (smrt_id %in% c('m64367e_221121_144411')){
                # Annamieke's brain
                print('#### This sample was Annamiekes brain sample')
                likely_match = tmp[1, ]
                likely_match[1, ] = NA
                likely_match$original_ID = 'Annamieks_sample'
            } else if (smrt_id %in% c('m64050_220124_203838', 'm64050_220123_094203', 'm64050_220225_130031', 'm64050_220226_235525', 'm64050_220228_090930', 'm64367e_220516_075357')){
                # Anke's brain
                print('#### This sample was Ankes brain sample')
                likely_match = tmp[1, ]
                likely_match[1, ] = NA
                likely_match$original_ID = 'Anke_sample'
            } else if (smrt_id %in% c('m64367e_220729_122853', 'm64050e_220729_122556', 'm64367e_220812_122636', 'm64050e_220812_122354')){
                # HG002
                print('#### This sample was HG002')
                likely_match = tmp[1, ]
                likely_match[1, ] = NA
                likely_match$original_ID = 'HG002'
            } else {
                # otherwise likely low coverage
                print('#### This sample was likely low coverage')
                likely_match = tmp[1, ]
                likely_match[1, ] = NA
            }
        } else if (nrow(likely_match) >1){
            print(paste0('!!! too many hits for ', f))
            # most times there's a duplication
            if (length(unique(likely_match$ID_100plus) == 1) || length(unique(likely_match$I_ID) == 1)){
                likely_match = likely_match[1, ]
                print('#### This was a duplicated sample.')
            }
        }
        likely_match$SMRT_ID = smrt_id
        all_checks_df = rbind(all_checks_df, likely_match)
    }

# 4. now we merge with coverage stats
    table(all_checks_df$SMRT_ID %in% coverage_stats$MOVIE_ID)
    all_check_df_coverage = merge(all_checks_df, coverage_stats, by.x = 'SMRT_ID', by.y = 'MOVIE_ID')

# 5. add project information and sequencing center
    all_check_df_coverage$SEQUENCING_CENTER = NA
    all_check_df_coverage$SEQUENCING_CENTER[grep('m64050', all_check_df_coverage$SMRT_ID)] = 'Amsterdam_Stitch'
    all_check_df_coverage$SEQUENCING_CENTER[grep('m64367', all_check_df_coverage$SMRT_ID)] = 'Amsterdam_Lilo'
    all_check_df_coverage$SEQUENCING_CENTER[is.na(all_check_df_coverage$SEQUENCING_CENTER)] = 'Radboud'
    # blood_brain_child project
    overview_genetics = fread('/project/holstegelab/Software/snakemake_pipeline/sample_check_data/20210615_overview_genetics.txt', h=T)
    bbc_project = overview_genetics[which(overview_genetics$PACBIO_SOMATIC == 'YES'),]
    all_check_df_coverage$PROJECT = ifelse(all_check_df_coverage$ID_GWAS %in% bbc_project$ID_GWAS, 'Blood_Brain_Child', NA)
    all_check_df_coverage$PROJECT[is.na(all_check_df_coverage$PROJECT)] = 'AD_Centenarians'
    all_check_df_coverage$PROJECT[which(all_check_df_coverage$diagnosis %in% c('Dementie_anders', 'NIID'))] = 'Other'
    all_check_df_coverage$PROJECT[is.na(all_check_df_coverage$diagnosis)] = 'Other'

# 6. save information
    cat('## Writing new freeze data..\n')
    all_check_df_coverage$PERC_HOMOLOGY = as.numeric(all_check_df_coverage$PERC_HOMOLOGY)
    outname = paste0('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/', Sys.Date(), '_freeze_sequencing_stats.txt')
    write.table(all_check_df_coverage, outname, quote=F, row.names=F, sep = "\t", dec = ',')