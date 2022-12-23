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
CCS_LAST="/home/holstegelab-ntesi/.conda/envs/smrtlink11/bin/ccs"
MD5="/usr/bin/md5sum"
EXTRACT_CCS="/project/holstegelab/Software/snakemake_pipeline/bin/extract_ccs_and_nonCCS.py"
EXTRACT_CCS_DEEPCONS="/project/holstegelab/Software/snakemake_pipeline/bin/extract_ccs_and_nonCCS_forDeepCons.py"
ALIGN="/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/pbmm2"
SAMPLE_CHECK="/project/holstegelab/Software/snakemake_pipeline/bin/sample_check.py"
MOSDEPTH="/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/mosdepth"
ASBT="/project/holstegelab/Software/nicco/tools/asbt/build/asbt"
ACTC='/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/actc'
DEEPCONSENSUS='/project/holstegelab/Software/conda/miniconda3_v1/envs/py37/bin/deepconsensus'
EXTRACT_READS_INFO='/project/holstegelab/Software/snakemake_pipeline/bin/extract_read_information.py'
PLOT_READS_INFO='/project/holstegelab/Software/snakemake_pipeline/bin/plot_read_information.R'

### RESOURCE PATHS
H38CCS='/project/holstegelab/Share/pacbio/resources/h38_ccs.mmi'
HG38SUBREADS='/project/holstegelab/Share/pacbio/resources/h38_subread.mmi'
CHM13CCS='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_hifi.mmi'
CHM13SUBREADS='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_subreads.mmi'
DCACHE_CONFIG='/project/holstegelab/Data/dcache.conf'
DEEPCONSENSUS_MODEL='/project/holstegelab/Software/nicco/tools/deepConsensun_examples/checkpoint'

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
        expand("{out_dir}/{out_name}_ccsLast/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)

# Rule to run CCS analysis on the cluster
rule ccs:
    input:
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    output:
        expand("{out_dir}/{out_name}_ccsLast/{out_name}.ccs.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name),
        expand("{out_dir}/{out_name}_ccsLast/{out_name}.ccs.report.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name),
        expand("{out_dir}/{out_name}_ccsLast/{out_name}.ccs.log", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    run:
        # check md5check and decide whether go on with the pipeline (if checksum was OK) or stop
        #RUN = open('%s/%s/%s/snakemake_checksum.txt' %(config["OUT_DIR"], folder_name, smrt_id)).readlines()[0].rstrip()
        RUN = 'True'
        if RUN == 'True':
                shell("{CCS_LAST} --min-rq 0.80 {input} {output} --report-file {out2} --log-file {out3}".format(CCS_LAST=CCS_LAST, input=input, output=output[0], out2=output[1], out3=output[2]))
