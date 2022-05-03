### PACKAGES
import os
import csv
import sys
import re

### SOFTWARE PATHS
SAMTOOLS='/project/holstegelab/Software/nicco/tools/samtools-1.11/samtools'
PYTHON39='/home/holstegelab-ntesi/.conda/envs/cpg/bin/python3.9'
PILEUP='/project/holstegelab/Share/oscar/software/CpG-main/scripts/RefAlnBam-to-ModsBed-SAMTags-modified.py'
SH='/usr/bin/sh'
DEEPVARIANT='/project/holstegelab/Software/nicco/bin/deepvariant/run_deepvariant_generic_forPipeline.sh'

### RESOURCES PATHS
HG38FA='/project/holstegelab/Share/pacbio/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa'
CHM13FA='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0.fa.gz'
MODELDIR="/project/holstegelab/Share/oscar/software/CpG-main/models/pileup_calling_model"

CACHE=None

### MAIN -- run before the rules -- identify input files and do md5sum check of the files
inp_bam_prefix = config["IN_FILES"].split(',')
out_bam_prefix = config["OUT_FILE"]

print("## Input files prefix --> %s" %(inp_bam_prefix))
print("## Output files prefix --> %s" %(out_bam_prefix))

# Rule all to check output files
rule all:
    input:
        # 1. merge hifi
        expand("{out_prefix}.merged.hifi.hg38.bam.bai", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.chm13.bam.bai", out_prefix = out_bam_prefix),

        # 2. merge non-hifi
        expand("{out_prefix}.merged.nonhifi.hg38.bam.bai", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.nonhifi.chm13.bam.bai", out_prefix = out_bam_prefix),

        # 3. deepvariant
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.vcf.gz", out_prefix = out_bam_prefix),
        #expand("{out_prefix}.merged.hifi.deepvariant.chm13.vcf.gz", out_prefix = out_bam_prefix)

        # 4. pileup analysis
        expand("{out_prefix}.met.hg38.combined.denovo.bed", out_prefix = out_bam_prefix),
        expand("{out_prefix}.met.chm13.combined.denovo.bed", out_prefix = out_bam_prefix)


### RULES FOR HG38
# Rule to merge hifi data aligned to hg38
rule merge_hifi_hg38:
    input:
        [x + '.ccs.primrose.hifi.hg38.bam' for x in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.hg38.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# Rule to index hifi data aligned to hg38
rule index_hifi_hg38:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.hg38.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# Rule to merge hifi data aligned to hg38
rule merge_nonhifi_hg38:
    input:
        [x + '.ccs.primrose.nonhifi.hg38.bam' for x in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.hg38.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# Rule to index hifi data aligned to hg38
rule index_nonhifi_hg38:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.nonhifi.hg38.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.hg38.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# Rule to perform pileup analysis for hg38
rule pileup_analysis_hg38:
    input:
        expand("{out_prefix}.merged.hifi.hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.hg38.bam.bai", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.met.hg38.combined.denovo.bed", out_prefix = out_bam_prefix)
    params:
        pfx = expand("{out_prefix}.met.hg38", out_prefix = out_bam_prefix)
    shell: """
        {PYTHON39} {PILEUP} -b {input[0]} -f {HG38FA} -o {params.pfx} -t 50 -p model -d {MODELDIR} --skip_bw_conversion
        """

# Rule to run deepvariant for hg38
rule deepvariant_hg38:
    input:
        expand("{out_prefix}.merged.hifi.hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.hg38.bam.bai", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.vcf.gz", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.gvcf.gz", out_prefix = out_bam_prefix)
    shell: """
        {SH} {DEEPVARIANT} {input[0]} {output[0]} {output[1]} {HG38FA}
        """

### RULES FOR CHM13
# Rule to merge hifi data aligned to chm13
rule merge_hifi_chm13:
    input:
        [x + '.ccs.primrose.hifi.chm13.bam' for x in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.chm13.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# Rule to index hifi data aligned to hg38
rule index_hifi_chm13:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.chm13.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# Rule to merge hifi data aligned to chm13
rule merge_nonhifi_chm13:
    input:
        [x + '.ccs.primrose.nonhifi.chm13.bam' for x in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.chm13.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# Rule to index hifi data aligned to hg38
rule index_nonhifi_chm13:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.nonhifi.chm13.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.chm13.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# Rule to perform pileup analysis for chm13
rule pileup_analysis_chm13:
    input:
        expand("{out_prefix}.merged.hifi.chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.chm13.bam.bai", out_prefix = out_bam_prefix)        
    output:
        expand("{out_prefix}.met.chm13.combined.denovo.bed", out_prefix = out_bam_prefix)
    params:
        pfx = expand("{out_prefix}.met.chm13", out_prefix = out_bam_prefix)
    shell: """
        {PYTHON39} {PILEUP} -b {input[0]} -f {CHM13FA} -o {params.pfx} -t 50 -p model -d {MODELDIR} --skip_bw_conversion
        """

# Rule to run deepvariant for hg38
rule deepvariant_chm13:
    input:
        expand("{out_prefix}.merged.hifi.chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.chm13.bam.bai", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.merged.hifi.deepvariant.chm13.vcf.gz", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.deepvariant.chm13.gvcf.gz", out_prefix = out_bam_prefix)
    shell: """
        {SH} {DEEPVARIANT} {input[0]} {output[0]} {output[1]} {CHM13FA}
        """
