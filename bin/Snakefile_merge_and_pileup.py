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
HIFIASM='~/.conda/envs/cpg/bin/hifiasm'
GFATOOLS='~/.conda/envs/cpg/bin/gfatools'
ALIGN='~/.conda/envs/cpg/bin/pbmm2'
BUSCO='~/.conda/envs/cpg/bin/busco'

### RESOURCES PATHS
HG38FA='/project/holstegelab/Share/pacbio/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa'
CHM13FA='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0.fa.gz'
MODELDIR="/project/holstegelab/Share/oscar/software/CpG-main/models/pileup_calling_model"
H38CCS='/project/holstegelab/Share/pacbio/resources/h38_ccs.mmi'
CHM13CCS='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_hifi.mmi'

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

        # 3. pileup analysis
        expand("{out_prefix}.met.hg38.combined.denovo.bed", out_prefix = out_bam_prefix),
        expand("{out_prefix}.met.chm13.combined.denovo.bed", out_prefix = out_bam_prefix),

        # 4. deepvariant
        expand("{out_prefix}.merged.hifi.deepvariant.hg38.vcf.gz", out_prefix = out_bam_prefix),
        #expand("{out_prefix}.merged.hifi.deepvariant.chm13.vcf.gz", out_prefix = out_bam_prefix),

        # 5. assemblies
        # convert merged hifi bam to fasta for assembly
        expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix),
        # hifiasm phased (convert gfa to fasta directly and align to both references)
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        # hifiasm primary (convert gfa to fasta directly and align)
        expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.bam", out_prefix = out_bam_prefix),

        # 6. busco analyses
        expand("{out_prefix}_primary_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_alternate_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap1_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap2_BUSCO/logs/busco.log", out_prefix = out_bam_prefix)

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

# Rule to convert BAM file to FASTA file
rule bam2fasta:
    input:
        expand("{out_prefix}.merged.hifi.hg38.bam", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix)
    shell: """
        {SAMTOOLS} fasta -@ 10 {input[0]} > {output[0]}
        """

# Rule to run hifiasm -- phased output
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

# Rule to run hifiasm -- primary contig with alternate contigs
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

# Rule to run Flye assembly

# Rule to convert GFA to fasta for alignment (haplotypes, primary, alternate)
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

# Alignment to GRCh38
rule align_hifiasm_hg38:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.bam", out_prefix = out_bam_prefix)
    shell: """
        {ALIGN} align --preset CCS {H38CCS} --unmapped --sort {input[0]} {output[0]} --log-level=INFO
        {ALIGN} align --preset CCS {H38CCS} --unmapped --sort {input[1]} {output[1]} --log-level=INFO
        {ALIGN} align --preset CCS {H38CCS} --unmapped --sort {input[2]} {output[2]} --log-level=INFO
        {ALIGN} align --preset CCS {H38CCS} --unmapped --sort {input[3]} {output[3]} --log-level=INFO
        """

# Alignment to CHM13
rule align_hifiasm_chm13:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg.fa", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap1.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.bp.hap2.p_ctg_chm13.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.bam", out_prefix = out_bam_prefix)
    shell: """
        {ALIGN} align --preset CCS {CHM13CCS} --unmapped --sort {input[0]} {output[0]} --log-level=INFO
        {ALIGN} align --preset CCS {CHM13CCS} --unmapped --sort {input[1]} {output[1]} --log-level=INFO
        {ALIGN} align --preset CCS {CHM13CCS} --unmapped --sort {input[2]} {output[2]} --log-level=INFO
        {ALIGN} align --preset CCS {CHM13CCS} --unmapped --sort {input[3]} {output[3]} --log-level=INFO
        """

# BUSCO analysis on the assembly fasta -- to be completed the output part once I know the output name
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