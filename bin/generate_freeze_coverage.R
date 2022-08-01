# Script to generate freeze of coverage information
#
# Briefly, takes coverage stats, and merge with additional information ragarding sample ID, sample type, project and sequencing center

# Libraries
    library(data.table)
    library(stringr)

# Main
# 1. read input coverage stats
    cat('## Reading coverage statistics..\n')
    coverage_stats = fread('/project/holstegelab/Share/pacbio/data_processed/coverage_smrt_cells.txt', h=T, stringsAsFactors=F, sep='\t')

# 2. read all files from the different projects
    cat('## Matching sample types and sequencing centers..\n')
    bbc_project = system("find /project/holstegelab/Share/pacbio/data_processed/blood_brain_child -name 'm*ccs*log'", intern = T)
    ad_chc_project = system("find /project/holstegelab/Share/pacbio/data_processed/ad_centenarians -name 'm*ccs*log'", intern = T)
    nijmegen = system("find /project/holstegelab/Share/pacbio/data_processed/nijmegen -name 'm*ccs*log'", intern = T)
    anke = system("find /project/holstegelab/Share/pacbio/data_processed/Anke_samples -name 'm*ccs*log'", intern = T)

# 3. parse all files and add columns such as project (ad-chc / bbc), sample and source
    # BBC
    bbc_project_tmp = data.frame(str_split_fixed(bbc_project, '/', 9), stringsAsFactors = F)
    bbc_project_tmp$X1 = NULL; bbc_project_tmp$NAME = str_replace_all(bbc_project_tmp$X9, '.ccs.log', '')
    bbc_project_tmp$source = 'VUMC'; bbc_project_tmp$source[grep('e_', bbc_project_tmp$NAME)] = 'NIJMEGEN'; bbc_project_tmp$source[grep('m64367e', bbc_project_tmp$NAME)] = "VUMC"
    bbc_project_tmp$date = str_split_fixed(bbc_project_tmp$NAME, '_', 3)[, 2]
    bbc_df = data.frame(data_path = bbc_project, project = bbc_project_tmp$X7, movie_id = bbc_project_tmp$NAME, sample_id = bbc_project_tmp$X8, source = bbc_project_tmp$source, date = bbc_project_tmp$date)
    # AD-CHC
    ad_chc_project_tmp = data.frame(str_split_fixed(ad_chc_project, '/', 8), stringsAsFactors = F)
    ad_chc_project_tmp$X1 = NULL; ad_chc_project_tmp$NAME = str_replace_all(ad_chc_project_tmp$X8, '.ccs.log', '')
    ad_chc_project_tmp$source = 'VUMC'; ad_chc_project_tmp$source[grep('e_', ad_chc_project_tmp$NAME)] = 'NIJMEGEN'; ad_chc_project_tmp$source[grep('m64367e', ad_chc_project_tmp$NAME)] = "VUMC"
    ad_chc_project_tmp$date = str_split_fixed(ad_chc_project_tmp$NAME, '_', 3)[, 2]
    ad_chc_df = data.frame(data_path = ad_chc_project, project = ad_chc_project_tmp$X7, movie_id = ad_chc_project_tmp$NAME, sample_id = rep(NA, nrow(ad_chc_project_tmp)), source = ad_chc_project_tmp$source, date = ad_chc_project_tmp$date)
    # NIJMEGEN
    nijmegen_tmp = data.frame(str_split_fixed(nijmegen, '/', 8), stringsAsFactors = F)
    nijmegen_tmp$X1 = NULL; nijmegen_tmp$NAME = str_replace_all(nijmegen_tmp$X8, '.ccs.log', '')
    nijmegen_tmp$source = 'VUMC'; nijmegen_tmp$source[grep('e_', nijmegen_tmp$NAME)] = 'NIJMEGEN'; nijmegen_tmp$source[grep('m64367e', nijmegen_tmp$NAME)] = "VUMC"; nijmegen_tmp$source[grep('m64050e', nijmegen_tmp$NAME)] = "VUMC"
    nijmegen_tmp$date = str_split_fixed(nijmegen_tmp$NAME, '_', 3)[, 2]
    nijmegen_df = data.frame(data_path = nijmegen, project = rep('ad_centenarians', nrow(nijmegen_tmp)), movie_id = nijmegen_tmp$NAME, sample_id = rep(NA, nrow(nijmegen_tmp)), source = nijmegen_tmp$source, date = nijmegen_tmp$date)
    # ANKE
    anke_tmp = data.frame(str_split_fixed(anke, '/', 8), stringsAsFactors = F)
    anke_tmp$X1 = NULL; anke_tmp$NAME = str_replace_all(anke_tmp$X8, '.ccs.log', '')
    anke_tmp$source = 'VUMC'; anke_tmp$source[grep('e_', anke_tmp$NAME)] = 'NIJMEGEN'; anke_tmp$source[grep('m64367e', anke_tmp$NAME)] = "VUMC"
    anke_tmp$date = str_split_fixed(anke_tmp$NAME, '_', 3)[, 2]
    anke_df = data.frame(data_path = anke, project = rep('anke_brain_bank', nrow(anke_tmp)), movie_id = anke_tmp$NAME, sample_id = rep(NA, nrow(anke_tmp)), source = anke_tmp$source, date = anke_tmp$date)

# 4. then add the coverage information
    all_runs = rbind(bbc_df, ad_chc_df, nijmegen_df, anke_df)
    all_runs_info = merge(all_runs, coverage_stats, by.x = 'movie_id', by.y = 'MOVIE_ID', all.x = T)

# 5. finally we need to add the actual sample match
    add_info = list()
    for (i in 1:nrow(all_runs_info)){
        # read sample match
        file_path = str_replace_all(all_runs_info$data_path[i], '.log', '.primrose.hifi.sample.txt')
        if (file.exists(file_path)){
            sample_info = fread(file_path, h=T, stringsAsFactors=F, sep="\t")
            tmp_info = sample_info[1, c('ID_GWAS', 'PERC_HOMOLOGY', 'SNPS_N', 'original_ID', 'ID_100plus', 'Study', 'diagnosis', 'sex', 'age')]
        } else {
            tmp_info = data.frame('ID_GWAS' = NA, 'PERC_HOMOLOGY' = NA, 'SNPS_N' = NA, 'original_ID' = NA, 'ID_100plus' = NA, 'Study' = NA, 'diagnosis' = NA, 'sex' = NA, 'age' = NA)
        }
        add_info[[(length(add_info) + 1)]] = tmp_info
    }
    add_info = rbindlist(add_info)
    # and add these information to the main dataframe
    all_runs_info = cbind(all_runs_info, add_info)
# 6. save information
    cat('## Writing new freeze data..\n')
    outname = paste0('/project/holstegelab/Share/pacbio/data_processed/coverage_freezes/', Sys.Date(), '_freeze_sequencing_stats.txt')
    write.table(all_runs_info, outname, quote=F, row.names=F, sep = "\t", dec = ',')