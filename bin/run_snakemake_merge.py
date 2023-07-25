# Script to guide merge pileup analysis on the samples from pacbio (non blood-brain-child project)
#
# The script prepares configuration files and launches merge_and_pileup analysis

# Libraries
import os
import pandas as pd
from os.path import exists
from datetime import date

# Functions
# function to read coverage information given a file path
def read_coverage(file_path):
    with open(file_path) as inpf:
        for line in inpf:
            if line.startswith('GLOBAL_MEDIAN_READ_LENGTH'):
                pass
            else:
                line = line.rstrip().split()
                smrt_id, global_cov = line[0], line[2]
    return smrt_id, global_cov

# function to read sample match information given a smrt cell id
def read_sample(file_path):
    fpath = file_path.replace('.ccs.primrose.hifi.hg38.coverage_summary.txt', '.ccs.primrose.hifi.sample.txt')
    sample_match = pd.read_csv(fpath, sep = '\t')
    # check if first is >0.94
    likely_match = sample_match.loc[sample_match['PERC_HOMOLOGY'] > 0.94]
    # if there's only 1 match with percentage of homology > 0.95, save the sample name and type
    # there are some exceptions for this: centenarian 100088 was duplicated but ok, anke's sample was not matched
    if likely_match.shape[0] == 1:
        id_gwas, diagnosis = likely_match.at[0, 'ID_GWAS'], likely_match.at[0, 'diagnosis']
    elif likely_match.shape[0] == 0:
        # put here exceptions solved
        if 'm64050e_230404_213248' in file_path:
            id_gwas, diagnosis = '18R0829', 'Probable_AD'
        else:
            print('!!! No matches for file %s' %(fpath))
            id_gwas, diagnosis = 'NA', 'NA'
    elif likely_match.at[0, 'ID_100plus'] == '100088':
        id_gwas, diagnosis = likely_match.at[0, 'ID_GWAS'], likely_match.at[0, 'diagnosis']        
    else:
        # put here exceptions solved
        if 'm64367e_230225_185330' in file_path:
            id_gwas, diagnosis = '9665575', 'Probable_AD'
        else:
            print('!!! Multiple matches for file %s' %(fpath))
            id_gwas, diagnosis = 'NA', 'NA'
    return id_gwas, diagnosis

# function to prepare configuration files for the merging operation
def prepare_config(smrt_cells, sample, diagnosis):
    if diagnosis == 'Centenarian':
        output_merged = '/project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/centenarian/' + sample
        tmp_dia = 'chc'
    elif diagnosis in ['Probable_AD', 'Possible_AD']:
        output_merged = '/project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/ad/' + sample
        tmp_dia = 'ad'
    else:
        output_merged = '/project/holstegelab/Share/pacbio/data_processed/merged_data_ad_centenarians/other/' + sample
        tmp_dia = 'other'
    configout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/config_%s_%s.yml' %(sample, tmp_dia), 'w')
    configout.write('IN_FILES : "%s"%s' %(','.join(smrt_cells), '\n\n'))
    configout.write('OUT_FILE : "%s"%s' %(output_merged, '\n'))
    configout.close()
    return None

# Main
# 1. list all files that finished the pipeline
all_files_ready_vu = os.popen('ls /project/holstegelab/Share/pacbio/data_processed/ad_centenarians/*coverage_summary*').read().rstrip().split('\n')
all_files_ready_radboud = os.popen('ls /project/holstegelab/Share/pacbio/data_processed/nijmegen/*coverage_summary*').read().rstrip().split('\n')
all_files_ready_other = os.popen('ls /project/holstegelab/Share/pacbio/data_processed/Anke_samples/*coverage_summary*').read().rstrip().split('\n')

# 2. combine all together
all_files = all_files_ready_vu + all_files_ready_radboud + all_files_ready_other

# 3. for each file, open the sample check, check if everything went fine and create a dictionary with sample ID and the relative smrt cells (path to data)
samples_files = {}
for f in all_files:
    # read coverage
    smrt_id, coverage = read_coverage(f)
    # read sample match
    try:
        id_gwas, diagnosis = read_sample(f)
    except:
        id_gwas, diagnosis = "NA", "NA"
    # save to dictionary
    if id_gwas != "NA":
        if id_gwas in samples_files.keys():
            samples_files[id_gwas].append([f, smrt_id, coverage, diagnosis])
        else:
            samples_files[id_gwas] = [[f, smrt_id, coverage, diagnosis]]

# 4. loop across samples and add the total coverage and number of smrt cells
for sample in samples_files.keys():
    smrt_cells_n = len(samples_files[sample])
    total_coverage = sum([float(x[-2]) for x in samples_files[sample]])
    samples_files[sample].append(smrt_cells_n)
    samples_files[sample].append(total_coverage)

# 5. for the samples with coverage >= 12 and at least 2 smrt cells, prepare configuration files
# also list all samples that have already been done
done_samples = [x.rstrip() for x in os.popen('ls /project/holstegelab/Software/snakemake_pipeline/config/config_merge/').readlines()]
done_samples = [x.split('_')[1] for x in done_samples]
to_do_samples = []
for sample in samples_files.keys():
    if sample != 'NA':
        if sample in done_samples:
            print('found %s' %(sample))
        else:
            smrt_cells_n, total_coverage = samples_files[sample][-2], samples_files[sample][-1]
            diagnosis = samples_files[sample][0][-1]
            if smrt_cells_n >= 2 and total_coverage >= 12:
                to_do_samples.append(sample)
                smrt_cells = [x[0].replace('.ccs.primrose.hifi.hg38.coverage_summary.txt', '') for x in samples_files[sample] if isinstance(x, list)]
                prepare_config(smrt_cells, sample, diagnosis)

# 6. finally submit jobs to the cluster (open a screen window, load conda environment, submit command)
# list all config files
config_files = os.popen("ls /project/holstegelab/Software/snakemake_pipeline/config/config_merge/ | egrep 'chc|ad'").read().rstrip().split('\n')
# to make sure we don't run the same sample multiple times, load file with submitted runs
submitted = os.popen('cat /project/holstegelab/Software/snakemake_pipeline/config/config_merge/merged_submitted.txt').read().rstrip().split()
new_submission = []
counter = 0
for config in config_files:
    # check if this submission was already done
    if config in submitted:
        print('!!! sample with config file --> %s was already submitted. Skipping to next sample.' %(config))
        pass
    else:
        counter = counter + 1
        if counter <15:
            print('XXX submitting sample with config file --> %s' %(config))
            # if the run is new, add it to the file
            fout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/merged_submitted.txt', 'a')
            fout.write('%s\n' %(config))
            fout.close()
            new_submission.append(config)
            # create an interactive screen session for the merging script
            os.system("screen -dmS 'merge_%s' /bin/bash -i" %(config.split('_')[1]))
            # then run in this screen session to load the right conda environment
            os.system("screen -S 'merge_%s' -X stuff 'conda activate cpg^M'" %(config.split('_')[1]))
            # then run command to go to the right directory
            os.system("screen -S 'merge_%s' -X stuff 'cd /project/holstegelab/Software/snakemake_pipeline/slurms_outputs^M'" %(config.split('_')[1]))
            # finally run the actual snakemake script
            os.system("screen -S 'merge_%s' -X stuff 'sh /project/holstegelab/Software/snakemake_pipeline/bin/submit_merge_and_pileup.sh /project/holstegelab/Software/snakemake_pipeline/config/config_merge/%s^M'" %(config.split('_')[1], config))

# 7. it's good to add merging information to the main excel i'm maintaining with run information. Generate a freeze tha can be easily copied into excel
if exists('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt'):
    fout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', 'a')
else:
    fout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', 'w')
    header = 'DATE\tSAMPLE\tSMRT_CELLS\tOUTPUT_IN\tSMRT_CELL_N\tCOMBINED_COVERAGE\tDIAGNOSIS\n'
    fout.write(header)
for config in new_submission:
    # get sample id to gather information
    sample_id = config.split('_')[1]
    sample_info = samples_files[sample_id]
    # get date
    date_today = date.today(); date_today = date_today.strftime("%d/%m/%Y")
    # sample is sample_id
    # get smrt cells id and number
    smrt_cells = [x[1] for x in sample_info if isinstance(x, list)]; smrt_cells_n = len(smrt_cells)
    # get output directory
    out_dir = os.popen('grep OUT_FILE /project/holstegelab/Software/snakemake_pipeline/config/config_merge/%s' %(config)).read().rstrip().split()[-1].replace('"', '')
    # get combined coverage
    combined_cov = sample_info[-1]
    # finally get diagnosis
    diagnosis = sample_info[0][-1]
    # add these info to the file
    fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(date_today, sample_id, ','.join(smrt_cells), out_dir, smrt_cells_n, str(combined_cov).replace('.', ','), diagnosis))
fout.close()
# to add: 9362652 -- 9274198

##########################################################################
##########################################################################
# separately, take the samples with a single smrt cell and enough coverage
single_smrt_enough = {}
for sample in samples_files.keys():
    if len(samples_files[sample]) == 3:
        if samples_files[sample][-1] >12:
            single_smrt_enough[sample] = [item for sublist in samples_files[sample] for item in (sublist if isinstance(sublist, list) else [sublist])]
# make a dataframe out of this
single_smrt_enough_df = pd.DataFrame.from_dict(single_smrt_enough, orient='index', columns=['data_path', 'smrt_id', 'coverage', 'diagnosis', 'n_smrt_cells', 'total_coverage'])
single_smrt_enough_df['id_gwas'] = single_smrt_enough.keys()

# these samples should be overlapped with the samples overview. read sample overview in
samples_overview = pd.read_csv('/project/holstegelab/Software/snakemake_pipeline/samples_overview/data_20230512/20230522_summary_smrtcells.txt', sep='\t')
# make subset of this dataframe to keep samples with enough coverage
samples_overview = samples_overview[samples_overview['sample'].isin(list(single_smrt_enough_df['id_gwas']))].copy()
single_samples_to_do = single_smrt_enough_df[single_smrt_enough_df['id_gwas'].isin(list(samples_overview['sample']))].copy()
# sort by coverage in descending order
samples_overview = samples_overview.sort_values('coverage', ascending=False)
single_samples_to_do = single_samples_to_do.sort_values('coverage', ascending=False)
single_samples_to_do_list = single_samples_to_do.values.tolist()
# also extract samples that are currently running
all_screens = [x.rstrip().split()[0].split('.')[-1] for x in list(os.popen('screen -ls'))]

# we run 4 centenarians and 4 ad
count_chc = 0; count_ad = 0
for s in single_samples_to_do_list:
    sample = s[-1]
    pheno = s[-4]
    # check the number of submissions
    if pheno == 'Probable_AD':
        run = True if count_ad <4 else False
        # compare with running processes
        if run == True:
            run = True if 'merge_%s' %(sample) not in all_screens else False
    elif pheno in ['Centenarian', 'Control_100plus']:
        pheno = 'Centenarian'
        run = True if count_chc <4 else False
        # compare with running processes
        if run == True:
            run = True if 'merge_%s' %(sample) not in all_screens else False
    # check if we have to run
    if run == True:
        if pheno == 'Centenarian':
            count_chc += 1
        else:
            count_ad += 1
        # identify data
        sample_data = single_smrt_enough[sample]
        # identify smrt cells
        smrt_cells = [sample_data[0].replace('.ccs.primrose.hifi.hg38.coverage_summary.txt', '') for x in samples_files[sample] if isinstance(x, list)]
        # prepare config file
        prepare_config(smrt_cells, sample, pheno)
        # finally submit jobs to the cluster (open a screen window, load conda environment, submit command)
        config_files = os.popen("ls /project/holstegelab/Software/snakemake_pipeline/config/config_merge/ | egrep 'chc|ad'").read().rstrip().split('\n')
        # to make sure we don't run the same sample multiple times, load file with submitted runs
        submitted = os.popen('cat /project/holstegelab/Software/snakemake_pipeline/config/config_merge/merged_submitted.txt').read().rstrip().split()
        new_submission = []
        counter = 0
        for config in config_files:
            # check if this submission was already done
            if config in submitted:
                #print('!!! sample with config file --> %s was already submitted. Skipping to next sample.' %(config))
                pass
            else:
                counter = counter + 1
                if counter <15:
                    print('XXX submitting sample with config file --> %s' %(config))
                    # if the run is new, add it to the file
                    fout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/merged_submitted.txt', 'a')
                    fout.write('%s\n' %(config))
                    fout.close()
                    new_submission.append(config)
                    # create an interactive screen session for the merging script
                    os.system("screen -dmS 'merge_%s' /bin/bash -i" %(config.split('_')[1]))
                    # then run in this screen session to load the right conda environment
                    os.system("screen -S 'merge_%s' -X stuff 'conda activate cpg^M'" %(config.split('_')[1]))
                    # then run command to go to the right directory
                    os.system("screen -S 'merge_%s' -X stuff 'cd /project/holstegelab/Software/snakemake_pipeline/slurms_outputs^M'" %(config.split('_')[1]))
                    # finally run the actual snakemake script
                    os.system("screen -S 'merge_%s' -X stuff 'sh /project/holstegelab/Software/snakemake_pipeline/bin/submit_merge_and_pileup.sh /project/holstegelab/Software/snakemake_pipeline/config/config_merge/%s^M'" %(config.split('_')[1], config))

## To UPDATE below
# next time i run this, check also the code below to ensure we add the single-smrt cell samples to the list of merged
# 7. it's good to add merging information to the main excel i'm maintaining with run information. Generate a freeze tha can be easily copied into excel
if exists('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt'):
    fout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', 'a')
else:
    fout = open('/project/holstegelab/Software/snakemake_pipeline/config/config_merge/freeze_merged_submitted.txt', 'w')
    header = 'DATE\tSAMPLE\tSMRT_CELLS\tOUTPUT_IN\tSMRT_CELL_N\tCOMBINED_COVERAGE\tDIAGNOSIS\n'
    fout.write(header)
for config in new_submission:
    # get sample id to gather information
    sample_id = config.split('_')[1]
    sample_info = samples_files[sample_id]
    # get date
    date_today = date.today(); date_today = date_today.strftime("%d/%m/%Y")
    # sample is sample_id
    # get smrt cells id and number
    smrt_cells = [x[1] for x in sample_info if isinstance(x, list)]; smrt_cells_n = len(smrt_cells)
    # get output directory
    out_dir = os.popen('grep OUT_FILE /project/holstegelab/Software/snakemake_pipeline/config/config_merge/%s' %(config)).read().rstrip().split()[-1].replace('"', '')
    # get combined coverage
    combined_cov = sample_info[-1]
    # finally get diagnosis
    diagnosis = sample_info[0][-1]
    # add these info to the file
    fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(date_today, sample_id, ','.join(smrt_cells), out_dir, smrt_cells_n, str(combined_cov).replace('.', ','), diagnosis))
fout.close()

