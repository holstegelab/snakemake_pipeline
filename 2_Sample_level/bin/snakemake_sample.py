###################################################################################################
###################################################################################################
### OVERVIEW OF PIPELINE AT THE LEVEL OF THE SAMPLE
# INPUT IS A CONFIGURATION FILE (CONFIG_FILE.YML, GIVEN IN THE FOLDER) WITH:
# 1. PATH TO (COMMA-SEPARATED) LIST OF PROCESSED (WITH THE STEP1 PIPELINE) FILES. ONLY SMRT CELL ID SHOULD BE GIVEN
# 2. OUTPUT DIRECTORY

# THE PIPELINE WILL SUBMIT JOBS TO SLURM. RESOURCES ARE DECLARED IN THE CLUSTER_CONFIG_FILE.YML (PRESENT IN THE FOLDER).

# TO RUN THE SCRIPT, THE COMMAND IS:
# snakemake -s /PATH/TO/THIS_FILE.PY --latency-wait 60 --configfile /PATH/TO/CONFIG_FILE.YML --cluster "sbatch --ntasks {cluster.ntasks} -c {cluster.ncpupertask} --time {cluster.time}" --cluster-config /PATH/TO/CLUSTER_CONFIG.YML --jobs 1

# THE JOBS ARE:
# 1. COMBINE ALIGNED HIFI AND (SEPARATELY) NON-HIFI DATA FROM ALL DEFINED SMRT CELLS
# 2. DE NOVO ASSEMBLY WITH HIFIASM - PHASED
# 3. DE NOVO ASSEMBLY WITH HIFIASM - PRIMARY CONTIG
# 4. DE NOVO ASSEMBLY WITH FLYE
# 5. ALIGNMENT OF ASSEMBLED CONTIGS (HIFIASM, FLYE) TO HG38 AND CHM13
# 6. SV CALLING WITH PBSV
# 7. SV CALLING WITH SNIFFLES
# 8. VARIANT CALLING WITH DEEPVARIANT
# 9. BUSCO ANALYSIS OF COMPLETENESS (CURRENTLY DISABLED)
# 10. SINGLE SMRT CELLS UPLOAD ON DCACHE, AND MOVE SINGLE SMRT CELLS TO TRASHBIN FOLDER
###################################################################################################
###################################################################################################

### PACKAGES
import os
import csv
import sys
import re
import random

### CONDA ENVIRONMENTS
# cpg
# install conda environments using the yml files

### SOFTWARE PATHS -- ADAPT THESE TO YOUR SYSTEM
SAMTOOLS='path/to/samtools-1.11/samtools'
PYTHON39='path/to/python3.9'
PILEUP='path/to/CpG-main/scripts/RefAlnBam-to-ModsBed-SAMTags-modified.py'
SH='/usr/bin/sh'
DEEPVARIANT='path/to/deepvariant/run_deepvariant_generic_forPipeline.sh'
HIFIASM='path/to/hifiasm'
GFATOOLS='path/to/gfatools'
ALIGN='path/to/pbmm2'
BUSCO='path/to/busco'
FLYE='path/to/flye'
PBSV='path/to/pbsv'
SNIFFLES='path/to/sniffles'
MINIMAP2='path/to/minimap2'

### RESOURCES PATHS -- ADAPT THESE TO YOUR SYSTEM
# HG38 GENOME IN FASTA FORMAT
HG38FA='path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa'
# CHM13 GENOME IN FASTA AND FASTA COMPRESSED
CHM13FAGZ='path/to/chm13v2.0.fa.gz'
CHM13FA='path/to/chm13v2.0.fa'
# MODEL FOR METHYLATION CALLING
MODELDIR="path/to/CpG-main/models/pileup_calling_model"
# DCACHE CONFIGURATION FILE 
DCACHE_CONFIG='path/to/dcache_processed.conf'
# TRASHBIN DIRECTORY -- AT THE END OF THE ANALYSIS, ALL DATA WILL BE PLACED THERE
TRASHBIN_PATH='path/to/trashbin'

CACHE=None

### MAIN -- IDENTIFY INPUT FILES
inp_bam_prefix = config["IN_FILES"].split(',')
# ADJUST IN CASE A DEMULTIPLEXED SAMPLE IS IN HERE
inp_bam_prefix = [inp.replace('.primrose.hifi.hg38.coverage_summary.txt', '') if '.demultiplexed' in inp else inp for inp in inp_bam_prefix]
out_bam_prefix = config["OUT_FILE"]
sample_type = out_bam_prefix.split('/')[-2]
dcache_folder = ""
# SET DIRECTORY TO COPY DATA TO DCACHE AT THE END OF THE ANALYSIS
if sample_type in ['ad', 'centenarian']:
    dcache_folder = 'ad_centenarians'
else:
    dcache_folder = 'other'
# OVERVIEW OF INPUT AND OUTPUT DATA
print("## Input files prefix --> %s" %(inp_bam_prefix))
print("## Output files prefix --> %s" %(out_bam_prefix))

# RULE TO DEFINE STEPS AND INPUT FILES
rule all:
    input:
        # 1. MERGE HIFI DATA
        expand("{out_prefix}.merged.hifi.hg38.bam.bai", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.chm13.bam.bai", out_prefix = out_bam_prefix),

        # 2. MERGE NON-HIFI DATA
        expand("{out_prefix}.merged.nonhifi.hg38.bam.bai", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.nonhifi.chm13.bam.bai", out_prefix = out_bam_prefix),

        # 3. PILEUP ANALYSIS FOR METHYLATION DATA
        expand("{out_prefix}.met.hg38.combined.denovo.bed", out_prefix = out_bam_prefix),
        expand("{out_prefix}.met.chm13.combined.denovo.bed", out_prefix = out_bam_prefix),

        # 4. DEEPVARIANT VARIANT CALLING ANALYSIS
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.vcf.gz", out_prefix = out_bam_prefix),
        #expand("{out_prefix}.merged.hifi.deepvariant.chm13.vcf.gz", out_prefix = out_bam_prefix),

        # 5. ASSEMBLIES
        # HIFIASM PHASED (CONVERT GFA TO FASTA AND ALIGN TO HG38 AND CHM13)
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        # HIFIASM PRIMARY (CONVERT GFA TO FASTA AND ALIGN TO HG38 AND CHM13)
        expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        # FLYE ASSEMBLY AND ALIGNMENT
        expand("{out_prefix}.hifi.flye/assembly_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/assembly_chm13.bam", out_prefix = out_bam_prefix),

        # 6. BUSCO ANALYSES OF COMPLETENESS OF GENOMES
        expand("{out_prefix}_primary_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_alternate_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap1_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap2_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/{sample_name}_BUSCO/logs/busco.log", out_prefix = out_bam_prefix, sample_name = out_bam_prefix.split('/')[-1]),

        # 7. STRUCTURAL VARIANT CALLING
        # PBSV
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.svcall.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.svcall.vcf'),
        # SNIFFLES
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.snf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.snf'),

        # 8. DATA CLEAN UP
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.files_copied.txt'),
        expand("{trashbin_path}/{fname}.files_copied.txt", trashbin_path = TRASHBIN_PATH, fname = out_bam_prefix.split('/')[-1])

### RULES FOR HG38
# MERGE HIFI DATA -- HG38
rule merge_hifi_hg38:
    input:
        [inp + '.ccs.primrose.hifi.hg38.bam' if '.demultiplexed.' not in inp else inp + '.primrose.hifi.hg38.bam' for inp in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.hg38.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# INDEX HIFI DATA -- HG38
rule index_hifi_hg38:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.hg38.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# MERGE NON-HIFI DATA -- HG38
rule merge_nonhifi_hg38:
    input:
        [inp + '.ccs.primrose.nonhifi.hg38.bam' if '.demultiplexed.' not in inp else inp + '.primrose.nonhifi.hg38.bam' for inp in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.hg38.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# INDEX NON-HIFI DATA -- HG38
rule index_nonhifi_hg38:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.nonhifi.hg38.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.hg38.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# PILEUP ANALYSIS -- HG38
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

# DEEPVARIANT ANALYSIS -- HG38
rule deepvariant_hg38:
    input:
        expand("{out_prefix}.merged.hifi.hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.hg38.bam.bai", out_prefix = out_bam_prefix),
        expand("{out_prefix}.met.hg38.combined.denovo.bed", out_prefix = out_bam_prefix),
        expand("{out_prefix}.met.chm13.combined.denovo.bed", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.vcf.gz", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.gvcf.gz", out_prefix = out_bam_prefix)
    shell: """
        {SH} {DEEPVARIANT} {input[0]} {output[0]} {output[1]} {HG38FA}
        """

### RULES FOR CHM13
# MERGE HIFI DATA -- CHM13
rule merge_hifi_chm13:
    input:
        [inp + '.ccs.primrose.hifi.chm13.bam' if '.demultiplexed.' not in inp else inp + '.primrose.hifi.chm13.bam' for inp in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.chm13.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# INDEX HIFI DATA -- CHM13
rule index_hifi_chm13:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.hifi.chm13.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# MERGE NON-HIFI DATA -- CHM13
rule merge_nonhifi_chm13:
    input:
        [inp + '.ccs.primrose.nonhifi.chm13.bam' if '.demultiplexed.' not in inp else inp + '.primrose.nonhifi.chm13.bam' for inp in inp_bam_prefix]
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.chm13.bam')
    shell: """
        {SAMTOOLS} merge {output[0]} {input} --threads 10
        """

# INDEX NON-HIFI DATA -- CHM13
rule index_nonhifi_chm13:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.nonhifi.chm13.bam')
    output:
        expand("{out_prefix}", out_prefix = out_bam_prefix + '.merged.nonhifi.chm13.bam.bai')
    shell: """
        {SAMTOOLS} index {input[0]} -@ 4
        """

# PILEUP ANALYSIS -- CHM13
rule pileup_analysis_chm13:
    input:
        expand("{out_prefix}.merged.hifi.chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.chm13.bam.bai", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.met.chm13.combined.denovo.bed", out_prefix = out_bam_prefix)
    params:
        pfx = expand("{out_prefix}.met.chm13", out_prefix = out_bam_prefix)
    shell: """
        {PYTHON39} {PILEUP} -b {input[0]} -f {CHM13FAGZ} -o {params.pfx} -t 50 -p model -d {MODELDIR} --skip_bw_conversion
        """

# DEEPVARIANT ANALYSIS -- CHM13
rule deepvariant_chm13:
    input:
        expand("{out_prefix}.merged.hifi.chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.chm13.bam.bai", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.merged.hifi.deepvariant.chm13.vcf.gz", out_prefix = out_bam_prefix),
        expand("{out_prefix}.merged.hifi.deepvariant.chm13.gvcf.gz", out_prefix = out_bam_prefix)
    shell: """
        {SH} {DEEPVARIANT} {input[0]} {output[0]} {output[1]} {CHM13FAGZ}
        """

# CONVERT BAM TO FASTA
rule bam2fasta:
    input:
        expand("{out_prefix}.merged.hifi.hg38.bam", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix))
    shell: """
        {SAMTOOLS} fasta -@ 10 {input[0]} > {output[0]}
        """

# HIFIASM -- PHASED OUTPUT
rule hifiasm_phased:
    input:
        expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.gfa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.gfa", out_prefix = out_bam_prefix)
    params:
        pfx = expand("{out_prefix}.hifi.hifiasm", out_prefix = out_bam_prefix)
    shell: """
        {HIFIASM} -o {params.pfx} -t 40 --n-hap 2 {input[0]}
    """

# HIFIASM -- PRIMARY
rule hifiasm_primary:
    input:
        expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.gfa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg.gfa", out_prefix = out_bam_prefix)
    params:
        pfx = expand("{out_prefix}.hifi.hifiasm", out_prefix = out_bam_prefix)
    shell: """
        {HIFIASM} -o {params.pfx} -t 40 --n-hap 2 {input[0]} --primary
        """

# FLYE
rule flye_assembly:
    input:
        expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix)
    params:
        pfx = expand("{out_prefix}.hifi.flye", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.flye/assembly.fasta", out_prefix = out_bam_prefix)
    shell: """
        {FLYE} -g 3.1g --threads 40 --out-dir {params.pfx} --pacbio-hifi {input[0]}
        """

# ALIGN FLYE ASSEMBLY -- HG38
rule align_flye_hg38:
    input:
        expand("{out_prefix}.hifi.flye/assembly.fasta", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.flye/assembly_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        """

# ALIGN FLYE ASSEMBLY -- CHM13
rule align_flye_chm13:
    input:
        expand("{out_prefix}.hifi.flye/assembly.fasta", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.flye/assembly_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        """

# BUSCO ANALYSIS ON FLYE ASSEMBLY
rule busco_flye:
    input:
        expand("{out_prefix}.hifi.flye/assembly.fasta", out_prefix = out_bam_prefix)
    params:
        pfx_path = expand("{path}", path = out_bam_prefix + '.hifi.flye'),
        pfx = expand("{sample_name}_BUSCO", sample_name = out_bam_prefix.split('/')[-1])
    output:
        expand("{out_prefix}.hifi.flye/{sample_name}_BUSCO/logs/busco.log", out_prefix = out_bam_prefix, sample_name = out_bam_prefix.split('/')[-1])
    shell: """
        {BUSCO} -f -i {input[0]} -o {params.pfx} -m genome -c 20 --out_path {params.pfx_path}
        """

# CONVERT GFA TO FASTA FOR ALIGNMENT (HAPLOTYPES, PRIMARY, ALTERNATE)
rule gfa2fasta:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.gfa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.gfa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.gfa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg.gfa", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    shell: """
        {GFATOOLS} gfa2fa {input[0]} > {output[0]}
        {GFATOOLS} gfa2fa {input[1]} > {output[1]}
        {GFATOOLS} gfa2fa {input[2]} > {output[2]}
        {GFATOOLS} gfa2fa {input[3]} > {output[3]}
        """

# ALIGN HIFIASM HAPLOTYPES -- HG38
rule align_hifiasm_hg38_haps:
    input:
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_hg38.sam", out_prefix = out_bam_prefix)),
        temp(expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[1]} > {output[1]}
        """

# ALIGN HIFIASM PRIMARY -- HG38
rule align_hifiasm_hg38_primary:
    input:
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        """

# ALIGN HIFIASM ALTERNATE -- HG38
rule align_hifiasm_hg38_alternate:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        """

# ALIGN HIFIASM HAPLOTYPES -- CHM13
rule align_hifiasm_chm13_haps:
    input:
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_chm13.sam", out_prefix = out_bam_prefix)),
        temp(expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[1]} > {output[1]}
        """

# ALIGN HIFIASM PRIMARY -- CHM13
rule align_hifiasm_chm13_primary:
    input:
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        """

# ALIGN HIFIASM ALTERNATE -- CHM13
rule align_hifiasm_chm13_alternate:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        """

# CONVERT SAM TO BAM
rule convertSAM_BAM:
    input:
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_hg38.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_hg38.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_chm13.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_chm13.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/assembly_hg38.sam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/assembly_chm13.sam", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/assembly_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/assembly_chm13.bam", out_prefix = out_bam_prefix)
    shell:"""
        {SAMTOOLS} view -bS -@ 3 {input[0]} > {output[0]}
        {SAMTOOLS} view -bS -@ 3 {input[1]} > {output[1]}
        {SAMTOOLS} view -bS -@ 3 {input[2]} > {output[2]}
        {SAMTOOLS} view -bS -@ 3 {input[3]} > {output[3]}
        {SAMTOOLS} view -bS -@ 3 {input[4]} > {output[4]}
        {SAMTOOLS} view -bS -@ 3 {input[5]} > {output[5]}
        {SAMTOOLS} view -bS -@ 3 {input[6]} > {output[6]}
        {SAMTOOLS} view -bS -@ 3 {input[7]} > {output[7]}
        {SAMTOOLS} view -bS -@ 3 {input[8]} > {output[8]}
        {SAMTOOLS} view -bS -@ 3 {input[9]} > {output[9]}
        """

# BUSCO ANALYSIS ON HIFIASM ASSEMBLIES
rule busco_analysis:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    params:
        pfx_path = expand("{path}", path = '/'.join(out_bam_prefix.split('/')[:-1])),
        pfx_p = expand("{sample_name}_primary_BUSCO", sample_name = out_bam_prefix.split('/')[-1]),
        pfx_a = expand("{sample_name}_alternate_BUSCO", sample_name = out_bam_prefix.split('/')[-1]),
        pfx_h1 = expand("{sample_name}_hap1_BUSCO", sample_name = out_bam_prefix.split('/')[-1]),
        pfx_h2 = expand("{sample_name}_hap2_BUSCO", sample_name = out_bam_prefix.split('/')[-1])
    output:
        expand("{out_prefix}_primary_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_alternate_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap1_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap2_BUSCO/logs/busco.log", out_prefix = out_bam_prefix)
    shell: """
        {BUSCO} -f -i {input[0]} -o {params.pfx_p} -m genome -c 20 --out_path {params.pfx_path}
        {BUSCO} -f -i {input[1]} -o {params.pfx_a} -m genome -c 20 --out_path {params.pfx_path}
        {BUSCO} -f -i {input[2]} -o {params.pfx_h1} -m genome -c 20 --out_path {params.pfx_path}
        {BUSCO} -f -i {input[3]} -o {params.pfx_h2} -m genome -c 20 --out_path {params.pfx_path}
        """

# SV CALLING STEP 1 -- FIND SIGNATURES -- PBSV
rule sv_discover_pbsv_reads:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.bam'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.bam.bai'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.bam'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.bam.bai')
    params:
        pfx_name = expand("{sample_prefix}", sample_prefix = out_bam_prefix.split('/')[-1])
    output:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.svsig.gz'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.svsig.gz'),
    shell: """
        {PBSV} discover --ccs -s {params.pfx_name} {input[0]} {output[0]}
        {PBSV} discover --ccs -s {params.pfx_name} {input[2]} {output[1]}
        """

# SV CALLING STEP 2 -- CALL SV -- PBSV
rule sv_call_pbsv_reads:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.svsig.gz'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.svsig.gz'),
    output:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.svcall.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.svcall.vcf'),
    shell: """
        {PBSV} call {HG38FA} --ccs -j 10 {input[0]} {output[0]}
        {PBSV} call {CHM13FA} --ccs -j 10 {input[1]} {output[1]}
        """

# SV CALLING -- SNIFFLES
rule sv_call_sniffles_reads:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.bam'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.bam.bai'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.bam'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.bam.bai')
    output:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.snf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.snf')
    shell: """
        {SNIFFLES} --input {input[0]} --vcf {output[0]} --snf {output[1]} --minsupport 2 --mapq 20 -t 10
        {SNIFFLES} --input {input[2]} --vcf {output[2]} --snf {output[3]} --minsupport 2 --mapq 20 -t 10
        """

# COPY SINGLE SMRT CELL DATA TO DCACHE ONCE ANALYSIS IS FINISHED
rule rclone_copy:
    input:
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.vcf.gz", out_prefix = out_bam_prefix),
        #expand("{out_prefix}_primary_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        #expand("{out_prefix}_alternate_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        #expand("{out_prefix}_hap1_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        #expand("{out_prefix}_hap2_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        #expand("{out_prefix}.hifi.flye/{sample_name}_BUSCO/logs/busco.log", out_prefix = out_bam_prefix, sample_name = out_bam_prefix.split('/')[-1]),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.svcall.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.svcall.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.vcf')
    output:
        #expand("{dcache_path}/{pheno_folder}/{fname}", dcache_path = DCACHE_PATH, pheno_folder = dcache_folder, fname = files_to_copy_names)
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.files_copied.txt')
    run:
        files_to_copy = ""
        for s in inp_bam_prefix:
            if '.demultiplexed.' in s:
                tmp_files = os.popen('ls ' + s + '.primrose.*').read().replace('\n', ',')
            else:
                tmp_files = os.popen('ls ' + s + '.ccs.*').read().replace('\n', ',')
            files_to_copy = files_to_copy + tmp_files
        files_to_copy = [x for x in files_to_copy.split(',') if x != '']
        for f in files_to_copy:
            # rclone copy
            shell('sleep 10s')
            shell('rclone copy -vv --progress --multi-thread-streams 1 --config {DCACHE_CONFIG} {input_file} dcache_processed:ccs/{dcache_folder_path}/'.format(input_file=f, DCACHE_CONFIG=DCACHE_CONFIG, dcache_folder_path=dcache_folder))
            shell('echo {input_file} >> {output_name}'.format(input_file=f, output_name=output))

# CHECK IF COPY TO DCACHE WENT OK AND MOVE DATA TO THE TRASHBIN 
rule rclone_check_move:
    input:
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.files_copied.txt')
    output:
        expand("{trashbin_path}/{fname}.files_copied.txt", trashbin_path = TRASHBIN_PATH, fname = out_bam_prefix.split('/')[-1])
    run:
        # read file containing copied files
        copied_files = open(out_bam_prefix + '.files_copied.txt').readlines()
        for f in copied_files:
            f = f.rstrip()
            # rclone check
            tmp_check = 'tmp_check_' + str(random.randint(1, 1000000)) + '.txt'
            shell('rclone check --config {DCACHE_CONFIG} {input_file} dcache_processed:ccs/{dcache_folder_path}/ --combined {tmp_check_path} --size-only'.format(input_file=f, DCACHE_CONFIG=DCACHE_CONFIG, dcache_folder_path=dcache_folder, tmp_check_path=tmp_check))
            tmp_check_open = open(tmp_check).readlines()[0].split(' ')[0]
            if tmp_check_open == '=':
                shell('rm {tmp_check_file}'.format(tmp_check_file=tmp_check))
                shell('mv {input_file} /project/holstegelab/Software/snakemake_pipeline/trashbin/'.format(input_file=f))
        # at the end, also move the file containing the files to copy
        shell('cp {checkfile} /project/holstegelab/Software/snakemake_pipeline/trashbin/'.format(checkfile = input))
