### PACKAGES
import os
import csv
import sys
import re
from os import path

### SOFTWARE PATHS
PYTHON="/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/python"
PRIMROSE="/project/holstegelab/Software/conda/miniconda3_v1/envs/py39/bin/primrose"
CCS="/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/ccs"
MD5="/usr/bin/md5sum"
EXTRACT_CCS="/project/holstegelab/Software/snakemake_pipeline/bin/extract_ccs_and_nonCCS.py"
EXTRACT_CCS_DEEPCONS="/project/holstegelab/Software/snakemake_pipeline/bin/extract_ccs_and_nonCCS_forDeepCons.py"
ALIGN="/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/pbmm2"
SAMPLE_CHECK="/project/holstegelab/Software/snakemake_pipeline/bin/sample_check.py"
MOSDEPTH="/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/mosdepth"
ASBT="/project/holstegelab/Software/nicco/tools/asbt/build/asbt"
SAMTOOLS='/project/holstegelab/Software/nicco/tools/samtools-1.11/samtools'
ACTC='/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/actc'
DEEPCONSENSUS='/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/deepconsensus'
LIMA='/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/lima'

### RESOURCE PATHS
H38CCS='/project/holstegelab/Share/pacbio/resources/h38_ccs.mmi'
HG38SUBREADS='/project/holstegelab/Share/pacbio/resources/h38_subread.mmi'
CHM13CCS='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_hifi.mmi'
CHM13SUBREADS='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_subreads.mmi'
DCACHE_CONFIG='/project/holstegelab/Data/dcache.conf'
DEEPCONSENSUS_MODEL='/project/holstegelab/Software/nicco/tools/deepConsensun_examples/checkpoint'
BARCODES='/project/holstegelab/Share/pacbio/data_processed/HLA_project/barcodes_sets/SMRTbell_Barcoded_Adapter_Plate_3.0_bc2001-bc2096.fasta'

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

### MAIN -- run before the rules
print("## Input subreads in dcache --> %s" %(config["IN_DIR"]))
print("## Output directory --> %s" %(config["OUT_DIR"]))
# create folder to place files copied from dcache
folder_name, smrt_id, out_name = config["IN_DIR"].split('/')[-3::]
if not path.isdir('%s/%s' %(config["OUT_DIR"], folder_name)):
    os.system('mkdir %s/%s' %(config["OUT_DIR"], folder_name))
# also check if a folder with the smrt_id name is present
if not path.isdir('%s/%s/%s' %(config["OUT_DIR"], folder_name, smrt_id)):
    os.system('mkdir %s/%s/%s' %(config["OUT_DIR"], folder_name, smrt_id))
# then adjust the output name and the path to be used as input for the dcache copy
out_name = out_name.split('.')[0]
dcache_path = '/'.join(config["IN_DIR"].split('/')[0:-1])

# Rule to manage output files and order of snakemake processes
rule all:
    input:
        # 1. copy from dcache
        #expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name),
        # 2. md5 check sum
        #expand("{out_dir}/{folder_name}/{smrt_id}/snakemake_checksum.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id),
        # 3. do ccs analysis keeping all kinetics values
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 4. run primrose analysis -- commented out 
        #expand("{out_dir}/{out_name}.ccs.primrose.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        # 5. extract hifi and non-hifi reads
        #expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
        expand("{out_dir}/{out_name}.ccs.demultiplexed.lima.summary", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.demultiplexed.0--0.fq", out_dir = config["OUT_DIR"], out_name = out_name)

# Rule to copy data from dcache
rule dcache_copy:
    output:
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    params:
        inp = expand("%s" %(dcache_path)),
        pfx = expand("{out_dir}/{folder_name}/{smrt_id}", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id)
    shell: """
        rclone copy -vv --progress --multi-thread-streams 1 --config {DCACHE_CONFIG} dcache:{params.inp}/ {params.pfx}/
        """

# Rule to perform md5 sum check if this is available
rule md5_check:
    input:
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    output:
        expand("{out_dir}/{folder_name}/{smrt_id}/snakemake_checksum.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id)
    params:
        inp = expand("{out_dir}/{folder_name}/{smrt_id}", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id)
    run:
        RUN = checkSum("%s/%s/%s" %(config["OUT_DIR"], folder_name, smrt_id), MD5)

# Rule to run CCS analysis on the cluster
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
            #shell("{CCS} --min-passes 0 --min-rq 0 {input[1]} {output[0]} --report-file {output[1]} --log-file {output[2]} --log-level INFO --all-kinetics")
            shell("{CCS} {input[1]} {output[0]} --report-file {output[1]} --log-file {output[2]} --log-level INFO")

# Rule for methylation analysis with primrose
rule primrose:
    input:
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        temp(expand("{out_dir}/{out_name}.ccs.primrose.bam", out_dir = config["OUT_DIR"], out_name = out_name))
    shell: """
        {PRIMROSE} --min-passes 0 {input[0]} {output[0]} 
        """

# Rule for extracting CCS and non-CCS reads
rule extract_ccs_nonccs:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name),
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {PYTHON} {EXTRACT_CCS} {input[0]} {output[0]} {output[1]}
        """

# Rule for demultiplexing hifi only
rule demultiplexing:
    input:
        expand("{out_dir}/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], out_name = out_name)
        #expand("{out_dir}/{out_name}.ccs.primrose.hifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    params:
        pfx = expand("{out_dir}/{out_name}.ccs.demultiplexed.bam", out_dir = config["OUT_DIR"], out_name = out_name)
        #pfx = expand("{out_dir}/{out_name}.ccs.primrose.hifi.demultiplexed.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.demultiplexed.lima.summary", out_dir = config["OUT_DIR"], out_name = out_name)
        #expand("{out_dir}/{out_name}.ccs.primrose.hifi.demultiplexed.lima.summary", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {LIMA} --ccs --split --min-score 80 --min-end-score 50 --min-ref-span 0.75 --same --ignore-missing-adapters -j 10 {input[0]} {BARCODES} {params.pfx}
        """

# Rule for demultiplexing hifi only
rule demultiplexing_nonhifi:
    input:
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    params:
        pfx = expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.demultiplexed.bam", out_dir = config["OUT_DIR"], out_name = out_name)
    output:
        expand("{out_dir}/{out_name}.ccs.primrose.nonhifi.demultiplexed.lima.summary", out_dir = config["OUT_DIR"], out_name = out_name)
    shell: """
        {LIMA} --ccs --split --min-score 80 --min-end-score 50 --min-ref-span 0.75 --same --ignore-missing-adapters -j 10 {input[0]} {BARCODES} {params.pfx}
        """

# Rule for converting bam to fastq
rule convert_to_fastq:
    input:
        expand("{out_dir}/{out_name}.ccs.demultiplexed.lima.summary", out_dir = config["OUT_DIR"], out_name = out_name)
    params:
        path = expand("{out_dir}/", out_dir = config["OUT_DIR"])
    output:
        expand("{out_dir}/{out_name}.ccs.demultiplexed.0--0.fq", out_dir = config["OUT_DIR"], out_name = out_name)
    run:
        shell('ls {out_dir}/{out_name}.ccs.demultiplexed.*.bam > {out_dir}/{out_name}.filelist.txt'.format(out_dir = config["OUT_DIR"], out_name = out_name))
        flist = open(config["OUT_DIR"] + '/' + out_name + '.filelist.txt').readlines()
        for f in flist:
            outname = f.replace('.bam', '.fq')
            shell('{SAMTOOLS} bam2fq {inp_bam} > {out_fq}'.format(SAMTOOLS = SAMTOOLS, inp_bam = f, out_fq = outname))