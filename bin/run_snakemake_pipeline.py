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
# 1. list all files in dcache: both the raw subreads (dcache) and processed files (dcache_processed)
config_dcache = '/project/holstegelab/Data/dcache.conf'
config_dcache_processed = '/project/holstegelab/Data/dcache_processed.conf'
# read files to be processed and parse them
flist = [x.rstrip() for x in os.popen('rclone ls --config ' + config_dcache + ' dcache:/')]
flist_subreads_qc1 = [x for x in flist if 'subreads.bam' in x]
flist_subreads_qc2 = [x for x in flist_subreads_qc1 if '.pbi' not in x]
flist_subreads_qc3 = [x.split()[-1] for x in flist_subreads_qc2]

# 2. list all processed files
# ad-centenarians folder on disk
proces_data_ad_chc = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/pacbio/data_processed/ad_centenarians -name '*.ccs.log'")]
# nijmegen folder on disk
proces_data_nijmegen = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/pacbio/data_processed/nijmegen -name '*.ccs.log'")]
# anke's folder on disk
proces_data_anke = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/pacbio/data_processed/Anke_samples -name '*.ccs.log'")]
# other folder on dist
proces_data_other = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/pacbio/data_processed/other_samples -name '*.ccs.log'")]
# blood-brain-child folder on disk
proces_data_bbc = [x.rstrip() for x in os.popen("find /project/holstegelab/Share/pacbio/data_processed/blood_brain_child -name '*.ccs.log'")]
# processed data on dcache
proces_data_dcache_qc1 = [x.rstrip() for x in os.popen('rclone ls --config ' + config_dcache_processed + ' dcache_processed:ccs/')]
proces_data_dcache_qc2 = [x for x in proces_data_dcache_qc1 if 'ccs.log' in x]
proces_data_dcache_qc3 = [x.split('/')[-1] for x in proces_data_dcache_qc2]
all_proces_data = proces_data_ad_chc + proces_data_dcache_qc3 + proces_data_anke + proces_data_bbc + proces_data_nijmegen + proces_data_other

# 3. loop through dcache files and see what has been done and what needs to be done
done = []; to_be_done = []
for subread in flist_subreads_qc3:
    movie_id = subread.split('/')[-1].replace('.subreads.bam', '')
    if len(list(filter(lambda x:movie_id in x, all_proces_data))) >0:
        done.append([movie_id, subread])
    else:
        to_be_done.append([movie_id, subread])

# 4. gather info for the runs that need to be done -- change from here
for run in to_be_done:
    # do not grep GenDX runs
    if 'Gendx' not in run[-1] and 'HLA' not in run[-1] and not run[0].startswith('.'):
        size = os.path.getsize('/home/holstegelab-ntesi/dcache/' + run[-1])/1e9
        date = run[0].split('_')[1]
        run.append(size); run.append(date)

# put data into dataframe
df = pd.DataFrame(to_be_done, columns = ['movie_id', 'data_path', 'size_gb', 'date'])

# 5. which run should we not run? [80-100Gb] gives about 1x coverage of hifi. I would exclude everything smaller than 80Gb
df_filtered = df[df['size_gb'] >= 80]

# 6. then split old runs (before 2022, blood-brain-child project, that were not done from the other AD-chc project)
to_be_done_datapath = list(df_filtered['data_path'])
newer_to_do = list(filter(lambda x: re.search(r'202[23]', x), to_be_done_datapath))
df_filtered_newer = df_filtered[df_filtered['data_path'].isin(newer_to_do)]
# also order by date
df_filtered_newer = df_filtered_newer.sort_values(by=['date'])

# 7. convert dataframe to list of lists
smrt_to_run = df_filtered_newer.values.tolist()

# 8. now it's time to create the config files for these smrt cells
for smrt in smrt_to_run:
    config_fname = prepare_config(smrt[1])
    smrt.append(config_fname)

# 9. to make sure we don't run things twice, let's read also the currently running screen processes
all_screen = [x.rstrip().split('\t') for x in list(os.popen('screen -ls'))]
all_screen = [x[1].split('.')[1] for x in all_screen if len(x) >1]

# 10. finally submit the snakemake pipeline (bash script, 1 argument which is the configuration file)
#for smrt in smrt_to_run:
    smrt = smrt_to_run[9]
    print('XXX submitting sample with config file --> %s' %(smrt[0]))
    # create an interactive screen session for the merging script
    screen_name = smrt[-1].split('/')[-1].replace('.yml', '').replace('config_', '')
    if screen_name not in all_screen:
        print(screen_name)
        os.system("screen -dmS '%s' /bin/bash -i" %(screen_name))
        # then run in this screen session to load the right conda environment
        os.system("screen -S '%s' -X stuff 'conda activate py37^M'" %(screen_name))
        # then run command to go to the right directory
        os.system("screen -S '%s' -X stuff 'cd /project/holstegelab/Software/snakemake_pipeline/slurms_outputs^M'" %(screen_name))
        # finally run the actual snakemake script
        os.system("screen -S '%s' -X stuff 'sh /project/holstegelab/Software/snakemake_pipeline/bin/submit_copy_and_pipeline.sh %s^M'" %(screen_name, smrt[-1]))
    else:
        pass

# LOGS (OLD)
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
    # RUN STARTED FOR r64050e_20220909_114615 (4 smrt cells)
    # RUN STARTED FOR r64167e_20220808_135837 (4 smrt cells)  
    # RUN STARTED FOR r64167e_20220715_20220715 (4 smrt cells)
    # RUN STARTED FOR r64050e_20220916_094122 (4 smrt cells)
    # RUN STARTED FOR r64367e_20220926_134945 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221004_153104 (4 smrt cell)
    # RUN STARTED FOR r64037e_20220812_121211 (4 smrt cells)
    # RUN STARTED FOR r64367e_20220912_104546 (4 smrt cells)
    # RUN STARTED FOR r64367e_20220920_083142 (2 smrt cells, 2 very low output)
    # RUN STARTED FOR r64050e_20220923_090617 (2 smrt cells, 2 very low output)
    # RUN STARTED FOR r64050e_20220930_085949 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221007_124347 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221010_142232 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221017_105039 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221017_103521 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221024_101008 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221024_101735 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221101_093432 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221101_093855 (4 smrt cells)
    # RUN STARTED FOR r64037e_20220822_121303 (3 smrt cells)*
    # RUN STARTED FOR r64037e_20221027_131940 (1 smrt cells)***
    # RUN STARTED FOR r64346e_20220825_091344 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221108_120220 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20221108_120723 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221115_131805 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20221115_132247 (3 smrt cells)*
    # RUN STARTED FOR r64050e_20221121_135634 (4 smrt cells)
    # RUN STARTED FOR r64367e_20221121_143147 (5 smrt cells)
    # RUN STARTED FOR r64050e_20221128_122837 (3 smrt cells)*
    # RUN STARTED FOR r64037e_20220829_123419 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20221129_094857 (1 smrt cells)***
    # RUN STARTED FOR r64367e_20221202_131637 (2 smrt cells)**
    # RUN STARTED FOR r64050e_20221205_111011 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20221205_140936 (4 smrt cells)
    # RUN STARTED FOR r64050e_20221212_145516 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20221212_150136 (3 smrt cells)*
    # RUN STARTED FOR r64102e_20220902_123837 (4 smrt cells)
    # RUN STARTED FOR r64346e_20220902_124909 (2 smrt cells)**
    # RUN STARTED FOR r64102e_20221111_142430 (4 smrt cells)
    # RUN STARTED FOR r64050e_20230109_142526 (4 smrt cells)
    # RUN STARTED FOR r64346e_20220907_132300 (1 smrt cells)***
    # RUN STARTED FOR r64050e_20230130_132350 (4 smrt cells)
    # RUN STARTED FOR r64037e_20220906_083239 (2 smrt cells)**
    # RUN STARTED FOR r64367e_20230206_135236 (4 smrt cells)
    # RUN STARTED FOR r64167e_20220907_133156 (4 smrt cells)
    # RUN STARTED FOR r64037e_20220929_082247 (3 smrt cells)*
    # RUN STARTED FOR r64050e_20221219_120452 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20230123_144644 (4 smrt cells)
    # RUN STARTED FOR r64037e_20220919_103932 (4 smrt cells)
    # RUN STARTED FOR r64167e_20220919_105929 (4 smrt cells)
    # RUN STARTED FOR r64347e_20220930_105833 (4 smrt cells)
    # RUN STARTED FOR r64367e_20230221_112935 (4 smrt cells)
    # RUN STARTED FOR r64346e_20220930_105730 (4 smrt cells)
    # RUN STARTED FOR r64050e_20230306_121710 (4 smrt cells)
    # RUN STARTED FOR r64167e_20221014_122045 (4 smrt cells)
    # RUN STARTED FOR r64367e_20230103_133345 (3 smrt cells)*
    # RUN STARTED FOR r64367e_20230109_143252 (4 smrt cells)
    # RUN STARTED FOR r64367e_20230116_131445 (1 smrt cells)***
    # RUN STARTED FOR r64367e_20230214_133248 (4 smrt cells)
    # RUN STARTED FOR r64050e_20230214_132436 (4 smrt cells)
    # RUN STARTED FOR r64050e_20230221_123319 (4 smrt cells)
    # RUN STARTED FOR r64367e_20230308_135456 (4 smrt cells)
    # RUN STARTED FOR r64050e_20230314_124647 (4 smrt cells)
    # RUN STARTED FOR r64367e_20230314_125931 (4 smrt cells)
    # RUN STARTED FOR r64050e_20230321_141438 (2 smrt cells)**
    # RUN STARTED FOR r64367e_20230321_142350 (1 smrt cells)***
