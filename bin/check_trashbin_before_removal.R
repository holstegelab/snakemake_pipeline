###########################################
# SCRIPT TO CHECK DATA IN TRASHBIN FOLDER #
# BEFORE DATA GETS REMOVED.               #
#                                         #
# THE IDEA IS TO CHECK WHETHER DATA WAS   #
# COPIED TO DCACHE, CHECK WHETHER THE     #
# COPY WAS SUCCESSFUL, AND IN CASE        #
# REMOVE DATA.                            #
###########################################

# set working directory to the trashbin folder
setwd('/project/holstegelab/Software/snakemake_pipeline/trashbin/')
# list all files in trashbin folder
trashbin = system('ls /project/holstegelab/Software/snakemake_pipeline/trashbin/', intern=T)
# exclude 2 files: files_with_issues.txt and rclone_checks.txt
# these files are not in dcache and therefore should not be checked
trashbin = trashbin[which(trashbin != 'files_with_issues.txt')]
trashbin = trashbin[which(trashbin != 'rclone_checks.txt')]
# set the path to the configuration file
dcache_config = '/project/holstegelab/Data/dcache_processed.conf'

# iterate through files in trashbin folder
for (f in trashbin){
    # rclone copy
    system(paste0('rclone copy -vv --progress --multi-thread-streams 1 --config ', dcache_config, ' ', f, ' dcache_processed:ccs/ad_centenarians'))
    # rclone check
    system(paste0('rclone check --size-only --combined rclone_checks.txt --config ', dcache_config, ' ', f, ' dcache_processed:ccs/ad_centenarians'))
    # read rclone check 
    rclone_check = read.table('rclone_checks.txt', h=F, stringsAsFactors=F)
    if (rclone_check$V2 == f & rclone_check$V1 == "="){
        # we can remove the file here
        system(paste0('rm ', f))
    } else {
        system(paste0('echo "', f, '" >> files_with_issues.txt'))
    }
}

################# THE END #################

###########################################
# DUE TO THE MACAROON EXPIRATION LAST     #
# JANUARY, WE NEED TO CHECK ALSO OTHER    #
# FOLDERS TO MAKE SURE THERE'S NO DATA    #
# LOSS.                                   #
###########################################

# list the directories to look (this can be re-used in future if the same thing happens)
other_directories = c('/home/holstegelab-ntesi/old_dcache', '/home/holstegelab-ntesi/old_dcache_cache', '/home/holstegelab-ntesi/to_upload_dcache')
# iterate through the directories
for (dir in other_directories){
    # set working directory in the directory
    setwd(dir)
    # find all files in the directories
    all_files = system(paste0('find ', dir, ' -type f'), intern = T)
    # exclude vfsmeta as they are just a report
    if (dir == '/home/holstegelab-ntesi/old_dcache_cache'){ all_files = all_files[grep('/vfs/', all_files)] }
    # iterate over the files
    for (f in all_files){
        # rclone copy
        system(paste0('rclone copy -vv --progress --multi-thread-streams 1 --config ', dcache_config, ' ', f, ' dcache_processed:ccs/ad_centenarians'))
        # rclone check
        system(paste0('rclone check --size-only --combined rclone_checks.txt --config ', dcache_config, ' ', f, ' dcache_processed:ccs/ad_centenarians'))
        # read rclone check 
        rclone_check = read.table('rclone_checks.txt', h=F, stringsAsFactors=F)
        fname = strsplit(f, '/')[[1]][length(strsplit(f, '/')[[1]])]
        if (rclone_check$V2 == fname & rclone_check$V1 == "="){
            # we can remove the file here
            system(paste0('rm ', f))
        } else {
            system(paste0('echo "', f, '" >> files_with_issues.txt'))
        }
    }
}