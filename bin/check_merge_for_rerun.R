# SCRIPT TO CHECK WHETHER THE MERGE ALGORITHM NEEDS TO BE RE-RAN OR NOT
# BECAUSE SOME SAMPLES MAY BE ADDED LATER ON AND THEREFORE THERE WOULD
# BE NEEDED TO RE-RUN THE PIPELINE
########################################################################

# LIBRARIES
    library(data.table)
    library(stringr)

# FUNCTIONS

# MAIN
    # 1. read freeze containing the merging information
    freeze = fread('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', h=T, stringsAsFactors=F)
    # 2. main loop across the freeze entries
    freeze$CHECK_INPUTS = NA
    for (i in 1:nrow(freeze)){
        cat(paste0('** processing file ', i, '              \r'))
        # get sample id and phenotype
        sample_id = freeze$SAMPLE[i]
        sample_pheno = ifelse(freeze$DIAGNOSIS[i] == 'Probable_AD', 'ad', 'chc')
        # get input files from freeze
        freeze_inputs = unlist(strsplit(freeze$SMRT_CELLS[i], ','))
        # get configuration file
        config = read.table(paste0('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/config_', sample_id, '_', sample_pheno, '.yml'), h=F, stringsAsFactors=F)
        # split input files
        input_files = unlist(strsplit(config$V3[1], ','))
        input_files_correct = c()
        for (x in input_files){ input_files_correct = c(input_files_correct, unlist(strsplit(x, '/'))[[length(unlist(strsplit(x, '/')))]]) }
        # now do the check
        if (setequal(freeze_inputs, input_files_correct) == TRUE){
            freeze$CHECK_INPUTS[i] = "same inputs"
        } else if (length(input_files_correct) == length(freeze_inputs)){
            freeze$CHECK_INPUTS[i] = "same length but different inputs"
        } else if (length(input_files_correct) > length(freeze_inputs)){
            freeze$CHECK_INPUTS[i] = 'new input found, need to re-run'
        } else if (length(input_files_correct) < length(freeze_inputs)){
            freeze$CHECK_INPUTS[i] = 'input lost, need to re-run'
        }
    }
    table(freeze$CHECK_INPUTS)