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
# define how many shards to do
shard_n = 100
all_shards = [x + 1 for x in range(shard_n)]

# Rule to manage output files and order of snakemake processes
rule all:
    input:
        # 1. copy data to spider
        #expand("{out_dir}/{folder_name}/{smrt_id}/snakemake_checksum.txt", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id),
        # 2. do ccs analysis keeping all kinetics values
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards),
        # 3. align subreads to ccs
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.aln_subreads.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards),
        # 4. deep consensus
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.deepConsensus.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards)

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
        expand("{out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam", out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name)
    output:
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards),
    run:
        # check md5check and decide whether go on with the pipeline (if checksum was OK) or stop
        RUN = open('%s/%s/%s/snakemake_checksum.txt' %(config["OUT_DIR"], folder_name, smrt_id)).readlines()[0].rstrip()
        if RUN == 'True':
            for x in range(len(all_shards)):
                outname = output[x]
                shell("{CCS} --min-rq 0.88 --chunk={n}/{n_tot} {input} {outname}".format(CCS=CCS, n=all_shards[x], n_tot=shard_n, input=input[1], outname=outname))

# Rule to calculate coverage summary
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

# Rule to align the created set of CCS reads back to the subreads
rule align_ccs_subreads_deepcons:
    input:
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards)
    output:
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.aln_subreads.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards)
    run:
        for x in range(len(all_shards)):
            input_name = input[x]
            output_name = output[x]
            shell("{ACTC} -j 40 {out_dir}/{folder_name}/{smrt_id}/{out_name}.subreads.bam {input_ccs} {output_aln}".format(ACTC=ACTC, out_dir = config["OUT_DIR"], folder_name = folder_name, smrt_id = smrt_id, out_name = out_name, input_ccs=input_name, output_aln=output_name))

# Rule to run deepConsensus
rule deepconsensus_run:
    input:
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.aln_subreads.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards)
    output:
        expand("{out_dir}/deepconsensus/{out_name}.ccs_{shard}.deepConsensus.bam", out_dir = config["OUT_DIR"], out_name = out_name, shard = all_shards)
    run:
        for x in range(len(all_shards)):
            input_ccs = input[x].replace('.aln_subreads.bam', '.bam')
            input_aln = input[x]
            out_file = output[x]
            shell("{DEEPCONSENSUS} run --subreads_to_ccs={input_aln} --ccs_bam={input_ccs} --checkpoint={DEEPCONSENSUS_MODEL} --output={out_file}".format(DEEPCONSENSUS=DEEPCONSENSUS, input_aln=input_aln, input_ccs=input_ccs, DEEPCONSENSUS_MODEL=DEEPCONSENSUS_MODEL, out_file=out_file))
