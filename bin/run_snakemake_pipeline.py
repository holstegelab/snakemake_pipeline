# Script to list all files in dcache and check those that have been done.
#
# In addition, this prepare also configuration files for snakemake.
#########################################################################

# Libraries
import os
import re
import pandas as pd

# Functions
# function to prepare configuration files for the snakemake pipeline
def prepare_config(data_path):
    config_prefix = 'config_%s_%s.yml' %(data_path.split('/')[-3], data_path.split('/')[-2].split('_')[0])
    if 'nijmegen' in data_path:
        config_fname = '/project/holstegelab/Software/snakemake_pipeline/config/config_files_snakemake/config_files_nijmegen/' + config_prefix
        out_dir = '/project/holstegelab/Share/pacbio/data_processed/nijmegen'
    else:
        config_fname = '/project/holstegelab/Software/snakemake_pipeline/config/config_files_snakemake/config_files_ad_chc/' + config_prefix
        out_dir = '/project/holstegelab/Share/pacbio/data_processed/ad_centenarians'
    configout = open(config_fname, 'w')
    configout.write('IN_DIR : "%s"%s' %(data_path.replace('/home/holstegelab-ntesi/dcache/', ''), '\n\n'))
    configout.write('OUT_DIR : "%s"%s' %(out_dir, '\n'))
    configout.close()
    return config_fname

# Main
# 1. list all files in dcache
dcache_data = [x.rstrip() for x in os.popen("find ~/dcache/ -name '*subreads.bam'")]

# 2. list all processed files
proces_data = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/pacbio/data_processed -name '*.ccs.log'")]

# 3. loop through dcache files and see what has been done and what needs to be done
done = []; to_be_done = []
for subread in dcache_data:
    movie_id = subread.split('/')[-1].replace('.subreads.bam', '')
    if len(list(filter(lambda x:movie_id in x, proces_data))) >0:
        done.append([movie_id, subread])
    else:
        to_be_done.append([movie_id, subread])

# 4. gather info for the runs that need to be done
for run in to_be_done:
    size = os.path.getsize(run[-1])/1e9
    date = run[0].split('_')[1]
    run.append(size); run.append(date)

# put data into dataframe
df = pd.DataFrame(to_be_done, columns = ['movie_id', 'data_path', 'size_gb', 'date'])

# 5. which run should we not run? [80-100Gb] gives about 1x coverage of hifi. I would exclude everything smaller than 80Gb
df_filtered = df[df['size_gb'] >= 80]

# 6. then split old runs (before 2022, blood-brain-child project, that were not done from the other AD-chc project)
to_be_done_datapath = list(df_filtered['data_path'])
newer_to_do = list(filter(lambda x: re.search(r'2022', x), to_be_done_datapath))
df_filtered_newer = df_filtered[df_filtered['data_path'].isin(newer_to_do)]
# also order by date
df_filtered_newer = df_filtered_newer.sort_values(by=['date'])

# 7. convert dataframe to list of lists
smrt_to_run = df_filtered_newer.values.tolist()

# 8. now it's time to create the config files for these smrt cells
for smrt in smrt_to_run:
    config_fname = prepare_config(smrt[1])
    smrt.append(config_fname)

# 9. finally submit the snakemake pipeline (bash script, 1 argument which is the configuration file)
#for smrt in smrt_to_run:
    #
    # !!!!! BE CAREFUL WITH RUNNING THIS FOR LOOP!! BETTER TO DO IT 1 BY 1 TO NOT MESS THINGS UP
    # !!!!! THE PROBLEM IS THAT IF SOME DATA IS BEIGN COPIED TO DCACHE, IT SHOULD NOT BE PROCESSED
    #
    # RUN STARTED FOR r64050_20220422_120813 AND r64050e_20220527_133058 [0, 1, 2, 3, 4 of smrt_to_run
    # RUN STARTED FOR r64050e_20220603_134449 (all 4 smrt cells)
    # RUN STARTED FOR r64367e_20220603_134547 (all 4 smrt cells)
    # RUN STARTED FOR r64050e_20220624_101248 (all 4 smrt cells)
    # RUN STARTED FOR r64367e_20220624_101256 (all 4 smrt cells)
    # RUN STARTED FOR r64050e_20220704_072158 (1 smrt cell available)
    # RUN STARTED FOR r64367e_20220707_114315 (all 4 smrt cells)
    # RUN STARTED FOR r64367e_20220715_121520 (all 4 smrt cells)
    # RUN STARTED FOR r64367e_20220722_091625 (all 4 smrt cells)
    # RUN STARTED FOR r64367e_20220729_121556 (all 4 smrt cells)
    # RUN STARTED FOR r64367e_20220805_090626 (all 2 smrt cells)
    # RUN STARTED FOR r64037e_20220623_103005 (all 4 smrt cells) -- Nijmegen -- trying monkey patch
    # RUN STARTED FOR r64102e_20220627_135041 (all 4 smrt cells) -- Nijmegen -- monkey patch 
    # RUN STARTED FOR r64167e_20220706_115319 (all 4 smrt cells) -- Nijmegen
    # RUN STARTED FOR r64102e_20220714_112024 (all 4 smrt cells) -- Nijmegen
    # RUN STARTED FOR r64367e_20220812_115756 (all 4 smrt cells)
    # RUN STARTED FOR r64050e_20220812_115528 (all 4 smrt cells)
    # RUN STARTED FOR r64167e_20220721_095203 (3 smrt cells -- 1 failed) -- Nijmegen
    # RUN STARTED FOR r64050e_20220819_141927 (3 smrt cells -- 1 very low output) -- Nijmegen
    # RUN STARTED FOR r64367e_20220819_142207 (all 4 smrt cells)
    # RUN STARTED FOR r64037e_20220804_105221 (all 4 smrt cells) -- Nijmegen
    # RUN STARTED FOR r64367e_20220825_143439 (3 smrt cells -- 1 failed)
    # RUN STARTED FOR r64050e_20220902_124118 (4 smrt cells)
    # RUN STARTED FOR r64050e_20220909_114615 (2 smrt cells for now - 2 running)
    # RUN STARTED FOR r64167e_20220808_135837 (4 smrt cells, 1 still running)  
    # RUN STARTED FOR r64167e_20220715_20220715 (4 smrt cells)
    # RUN STARTED FOR r64050e_20220916_094122 (2 smrt cells)
    # RUN STARTED FOR r64367e_20220926_134945 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221004_153104 (1 smrt cell)
    # 
    #
    print('XXX submitting sample with config file --> %s' %(smrt[0]))
    # create an interactive screen session for the merging script
    screen_name = smrt[-1].split('/')[-1].replace('.yml', '').replace('config_', '')
    os.system("screen -dmS '%s' /bin/bash -i" %(screen_name))
    # then run in this screen session to load the right conda environment
    os.system("screen -S '%s' -X stuff 'conda activate py37^M'" %(screen_name))
    # then run command to go to the right directory
    os.system("screen -S '%s' -X stuff 'cd /project/holstegelab/Software/snakemake_pipeline/slurms_outputs^M'" %(screen_name))
    # finally run the actual snakemake script
    os.system("screen -S '%s' -X stuff 'sh /project/holstegelab/Software/snakemake_pipeline/bin/submit_copy_and_pipeline.sh %s^M'" %(screen_name, smrt[-1]))
