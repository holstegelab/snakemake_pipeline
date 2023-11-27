###################################################################################################
###################################################################################################
### OVERVIEW OF PIPELINE AT THE LEVEL OF THE SINGLE SMRT CELL
# INPUT IS A CONFIGURATION FILE (CONFIG_FILE.YML, GIVEN IN THE FOLDER) WITH:
# 1. INPUT SUBREADS DATA (FROM DCACHE)
# 2. OUTPUT DIRECTORY

# THE PIPELINE WILL SUBMIT JOBS TO SLURM. RESOURCES ARE DECLARED IN THE CLUSTER_CONFIG_FILE.YML (PRESENT IN THE FOLDER).

# TO RUN THE SCRIPT, THE COMMAND IS:
# snakemake -s /PATH/TO/THIS_FILE.PY --latency-wait 60 --configfile /PATH/TO/CONFIG_FILE.YML --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask} --time {cluster.time}" --cluster-config PATH/TO/CLUSTER_CONFIG.YML --jobs 1

# THE JOBS ARE:
# 1. DOWNLOAD FROM DCACHE
# 2. MD5 CHECKSUM
# 3. CCS ALGORITHM
# 4. SEPARATE HIFI FROM NON-HIFI DATA
# 5. ALIGNMENT (HIFI AND NON-HIFI) TO HG38 AND CHM13
# 6. SAMPLE CHECK (COMPARISON WITH GWAS DATA)
# 7. COVERAGE SUMMARY
# 8. READ LENGTH AND PASSES DISTRIBUTION PLOT
# 9. DEEPCONSENSUS (BETA VERSION)
###################################################################################################
###################################################################################################

### PACKAGES
import os
import csv
import sys
import re
from os import path

### CONDA ENVIRONMENTS
# py37 (python3.7 and most packages)
# py39 (python3.9 was added later as primrose specifically needed python3.9)
# install conda environments using the yml files
# THESE ENVIRONMENTS ARE GIVEN IN THE FOLDER, ARE LIKELY REDUNDANT AND ARE NOT BASED ON THE LATEST VERSION OF SOME PACKAGES (WE KEPT SOFTWARE VERSIONS THE SAME FOR COMPATIBILITY)

### SOFTWARE PATHS -- ADAPT THESE TO YOUR SYSTEM
PYTHON="usr/bin/python"
PRIMROSE="path/to/primrose"
CCS="path/to/ccs"
MD5="/usr/bin/md5sum"
ALIGN="path/to/pbmm2"
ASBT="path/to/asbt"
# ONLY IF YOU WILL USE DEEPCONSENSUS
DEEPCONSENSUS='path/to/deepconsensus'
ACTC='path/to/actc'

#### HOMEMADE SCRIPTS -- ADAPT PATHS IF NEEDED
# SCRIPTS TO SEPARATE HIFI FROM NON-HIFI READS
EXTRACT_CCS="/path/to/snakemake_pipeline/bin/extract_ccs_and_nonCCS.py"
EXTRACT_CCS_DEEPCONS="/path/to/snakemake_pipeline/bin/extract_ccs_and_nonCCS_forDeepCons.py"
# SCRIPT TO PERFORM SAMPLE CHECK COMPARED TO GWAS DATA
SAMPLE_CHECK="/path/to/snakemake_pipeline/bin/sample_check.py"
# SCRIPT TO EXTRACT READ INFORMATION (NUMBER OF PASSES AND READ QUALITY) AND PLOT
EXTRACT_READS_INFO='/path/to/snakemake_pipeline/bin/extract_read_information.py'
PLOT_READS_INFO='/path/to/snakemake_pipeline/bin/plot_read_information.R'

### RESOURCE PATHS
# INDEXED REFERENCE GENOME HG38 (HIFI AND NON-HIFI)
H38CCS='path/to/h38_ccs.mmi'
HG38SUBREADS='path/to/h38_subread.mmi'
# INDEXED REFERENCE GENOME CHM13 (HIFI AND NON-HIFI)
CHM13CCS='path/to/chm13v2.0_hifi.mmi'
CHM13SUBREADS='path/to/chm13v2.0_subreads.mmi'
# DCACHE CONFIGURATION FILE
DCACHE_CONFIG='path/to/dcache.conf'
# DEEPCONSENSUS MODEL
DEEPCONSENSUS_MODEL='path/to/deepConsensun_examples/checkpoint'

CACHE=None

### FUNCTIONS
# Function to perform MD5 checksum and return True/False (RUN)
def checkSum(path, MD5):
    # first check if any .md5 file is present in the input directory
    flist = os.popen("ls %s" %(path)).read().rstrip().split("\n")
    md5 = list(filter(lambda x: re.search(r'md5$', x), flist))
    match_done = list(filter(lambda x: re.search(r'md5checksum.txt', x), flist))
    if (len(md5) == 1) & (len(match_done) == 0):
        print("## MD5 file found. Now checking all files were correctly copied.")
        # read and store MD5 file
        md5_file = []
        with open(path + "/" + md5[0]) as inpmd5:
            for line in inpmd5:
                if line.startswith('export') or 'md5' in line:
                    pass
                else:
                    line = line.rstrip().split()
                    sumNumb, filename = line[0], line[1].split("/")[-1]
                    md5_file.append([sumNumb, filename])
        # write new output with full paths
        outf = open((path + "/" + "md5checksum.md5"), "w")
        for x in md5_file:
            outf.write("%s\t%s\n" %(x[0], (path + "/" + x[1])))
        outf.close()
        # then do the check
        os.system("%s -c %s > %s" %(MD5, (path + "/" + "md5checksum.md5"), (path + "/" + "md5checksum.txt")))
        # ok now let's check if everything went ok
        RUN = True
        with open((path + "/" + "md5checksum.txt")) as md5check:
            for line in md5check:
                if RUN == False:
                    pass
                else:
                    check = line.rstrip().split(" ")
                    if check[-1] != "OK":
                        RUN = False
        # finally delete the created (alias) md5 file
        os.system("rm %s" %((path + "/" + "md5checksum.md5")))
    elif len(match_done) == 1:
        print("## MD5 output found. Checking the file.")
        RUN = True
        with open((path + "/" + "md5checksum.txt")) as md5check:
            for line in md5check:
                if RUN == False:
                    pass
                else:
                    check = line.rstrip().split(" ")
                    if check[-1] != "OK":
                        RUN = False
    else:
        print("## MD5 file NOT found. Will assume all files were correctly copied.")
        RUN = True
    # finally create a file that will be used by snakemake to check outputs
    outf = open(path + '/snakemake_checksum.txt', 'w')
    outf.write('%s\n' %(RUN))
    outf.close()
    return RUN

### MAIN
print("## Input subreads in dcache --> %s" %(config["IN_DIR"]))
print("## Output directory --> %s" %(config["OUT_DIR"]))

# COPY THE SMRT CELL DATA TO LOCAL -- CREATE FOLDER WITH RUN ID
folder_name, smrt_id, out_name = config["IN_DIR"].split('/')[-3::]
if not path.isdir('%s/%s' %(config["OUT_DIR"], folder_name)):
    os.system('mkdir %s/%s' %(config["OUT_DIR"], folder_name))
# CHECK IF A FOLDER WITH THE SAME SMRT ID IS PRESENT
if not path.isdir('%s/%s/%s' %(config["OUT_DIR"], folder_name, smrt_id)):
    os.system('mkdir %s/%s/%s' %(config["OUT_DIR"], folder_name, smrt_id))
# ADJUST THE OUTPUT NAME AND PATHS TO BE USED FOR THE DCACHE COPY
out_name = out_name.split('.')[0]
dcache_path = '/'.join(config["IN_DIR"].split('/')[0:-1])

# RULE TO MANAGE WHICH STEPS WILL BE DONE (BASED ON OUTPUT FILE EXTENSION)
rule all:
    input:
        # 1. COPY FROM DCACHE TO LOCAL
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name),
        # 2. MD5 CHECK 
        expand("{out_dir}/{folder_name}/{smrt_id}/snakemake_checksum.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id),
        # 3. CCS ANALYSIS KEEPING ALL KINETICS INFORMATION
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 4. EXTRACT HIFI AND NON-HIFI DATA
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 5. ALIGN HIFI TO HG38
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.hg38.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 6. ALIGN NON-HIFI TO HG38
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.hg38.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 7. ALIGN HIFI TO CHM13
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.chm13.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 8. ALIGN NON-HIFI TO CHM13
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.chm13.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 9. SAMPLE CHECK
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.sample.txt", out_dir = config["OUT_DIR"], out_name = out_name),
        # 10. GENERATE COVERAGE SUMMARY AND PLOT
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.hg38.coverage_summary.txt", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.reads_summary.png", out_dir = config["OUT_DIR"], out_name = out_name),
        # 11. CREATE INPUT FOR DEEPCONSENSUS -- COMMENTED AS WE CURRENTLY DON'T DO THIS
        #expand("{out_dir}/{out_name}.ccs.hifi_for_deepCons.bam", out_dir = config["OUT_DIR"], out_name = out_name)
        # 12. RUN DEEPCONSENSUS -- COMMENTED AS WE CURRENTLY DON'T DO THIS
        #expand("{out_dir}/{out_name}.ccs.deepConsensus.bam", out_dir = config["OUT_DIR"], out_name = out_name)

# RULE TO DOWNLOAD DATA FROM DCACHE, STEPS ARE: STAGE -- DOWNLOAD -- UNSTAGE
rule dcache_copy:
    output:
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    params:
        inp = expand("%s" %(dcache_path)),
        pfx = expand("{out_dir}/{folder_name}/{smrt_id}", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id)
    shell: """
        ada --tokenfile {DCACHE_CONFIG} --stage {params.inp}
        sleep 7200
        rclone copy -vv --progress --multi-thread-streams 1 --config {DCACHE_CONFIG} dcache:{params.inp}/ {params.pfx}/
        ada --tokenfile {DCACHE_CONFIG} --unstage {params.inp}
        """

# RULE TO PERFORM MD5 CHECK IF THIS IS AVAILABLE
rule md5_check:
    input:
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    output:
        expand("{out_dir}/{folder_name}/{smrt_id}/snakemake_checksum.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id)
    params:
        inp = expand("{out_dir}/{folder_name}/{smrt_id}", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id)
    run:
        RUN = checkSum("%s/%s/%s" %(config["OUT_DIR"], folder_name, smrt_id), MD5)

# RULE TO PERFORM CCS ANALYSIS
rule ccs:
    input:
        expand("{out_dir}/{folder_name}/{smrt_id}/snakemake_checksum.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id),
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name),
    output:
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.report.txt", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.log", out_dir = config["OUT_DIR"], out_name = out_name)
    run:
        # check md5check and decide whether go on with the pipeline (if checksum was OK) or stop
        RUN = open('%s/%s/%s/snakemake_checksum.txt' %(config["OUT_DIR"], folder_name, smrt_id)).readlines()[0].rstrip()
        if RUN == 'True':
            shell("{CCS} --min-passes 0 --min-rq 0 {input[1]} {output[0]} --report-file {output[1]} --log-file {output[2]} --log-level INFO --all-kinetics")

# RULE FOR METHYLATION ANALYSIS WITH PRIMROSE
rule primrose:
    input:
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        temp(expand("{out_dir}/{out_name}.ccs.primrose.bam", out_dir = config["OUT_DIR"], out_name = out_name))
    shell: """
        {PRIMROSE} --min-passes 0 {input[0]} {output[0]} 
        """

# RULE TO EXTRACT HIFI AND NON-HIFI DATA
rule extract_ccs_nonccs:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {PYTHON} {EXTRACT_CCS} {input[0]} {output[0]} {output[1]}
        """

# RULE TO ALIGN HIFI DATA TO HG38
rule align_ccs:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.hg38.bam", out_dir = config["OUT_DIR"], out_name = out_name),
    shell: """
        {ALIGN} align --preset CCS {H38CCS} --unmapped --sort {input[0]} {output[0]} --log-level=INFO
        """

# RULE TO ALIGN NON-HIFI DATA TO HG38
rule align_nonccs:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.hg38.bam", out_dir = config["OUT_DIR"], out_name = out_name),
    shell: """
        {ALIGN} align --preset SUBREAD {HG38SUBREADS} --unmapped --sort {input[0]} {output[0]} --log-level=INFO
        """

# RULE TO ALIGN HIFI DATA TO CHM13
rule align_ccs_chm13:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.chm13.bam", out_dir = config["OUT_DIR"], out_name = out_name),
    shell: """
        {ALIGN} align --preset CCS {CHM13CCS} --unmapped --sort {input[0]} {output[0]} --log-level=INFO
        """

# RULE TO ALIGN NON-HIFI DATA TO CHM13
rule align_nonccs_chm13:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.chm13.bam", out_dir = config["OUT_DIR"], out_name = out_name),
    shell: """
        {ALIGN} align --preset SUBREAD {CHM13SUBREADS} --unmapped --sort {input[0]} {output[0]} --log-level=INFO
        """

# RULE TO DO SAMPLE CHECK
rule sample_check:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.hg38.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.sample.txt", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {PYTHON} {SAMPLE_CHECK} {input[0]} {output[0]}
        """

# RULE TO CALCULATE COVERAGE STATISTICS
rule coverage_summary:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.hg38.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.hg38.coverage_summary.txt", out_dir = config["OUT_DIR"], out_name = out_name)
    params:
        pfx_name = expand("{out_name}", out_name = out_name),
        pfx_main = expand("/project/holstegelab/Share/pacbio/data_processed/coverage_smrt_cells.txt")
    shell: """
        printf "GLOBAL_MEDIAN_READ_LENGTH\tGLOBAL_COVERAGE\tMAPPED_COVERAGE\tMAPPED_READS\tALT_COVERAGE\tALT_READS\tUNMAPPED_COVERAGE\tUNMAPPED_READS\n{params.pfx_name}\t" > {output[0]}
        {ASBT} cov -g 3088000000 {input[0]} >> {output[0]}
        tail -1 {output[0]} >> {params.pfx_main}
        """

# RULE TO PLOT READ-LENGTH AND PASSES DISTRIBUTION IN HIFI AND NON-HIFI DATA
rule read_length_distr:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.reads_summary.txt", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.reads_summary.png", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {PYTHON} {EXTRACT_READS_INFO} {input[0]} {input[1]}
        /usr/bin/Rscript {PLOT_READS_INFO} {output[0]}
        """

###################################################################################################
# DEEPCONSENSUS -- NOT PERFORMED AT THE MOMENT, BUT FUNCTIONAL
# DEEPCONSENSUS IS A 2-STEP PROCESS: FIRST, CCS SHOULD BE RUN WITH --min-rq=0.88 AND DEFAULT FOR PASSES.
# THEN, DEEPCONSENSUS FIRST ALIGNS SUBREADS TO THE ACTUAL CCS (ACTC), AND THEN DEEPCONSENSUS POLISHES THE READS
# THUS, WE NEED TO HAVE:
# 1. THE RAW SUBREADS DATA
# 2. THE HIFI DATA WITH REFINED PARAMETERS

# RULE TO EXTRACT HIFI READS FOR DEEPCONSENSUS
rule extract_ccs_deepcons:
    input:
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.hifi_for_deepCons.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {PYTHON} {EXTRACT_CCS_DEEPCONS} {input[0]} {output[0]}
        """

# RULE TO ALIGN SUBREADS TO THE NEW SET OF HIFI READS
rule align_ccs_subreads_deepcons:
    input:
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.hifi_for_deepCons.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.hifi_for_deepCons_aln_subreads.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {ACTC} -j 40 {input[0]} {input[1]} {output[0]}
        """

# RULE TO RUN DEEPCONSENSUS
rule deepconsensus_run:
    input:
        expand("{out_dir}/{out_name}.ccs.hifi_for_deepCons_aln_subreads.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.hifi_for_deepCons.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.deepConsensus.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {DEEPCONSENSUS} run --subreads_to_ccs={input[0]} --ccs_bam={input[1]} --checkpoint={DEEPCONSENSUS_MODEL} --output={output[0]}
        """
