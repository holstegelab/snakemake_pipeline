### PACKAGES
import os
import csv
import sys
import re
import random

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
FLYE='~/.conda/envs/cpg/bin/flye'
PBSV='~/.conda/envs/cpg/bin/pbsv'
SNIFFLES='~/.conda/envs/cpg/bin/sniffles'
MINIMAP2='~/.conda/envs/cpg/bin/minimap2'

### RESOURCES PATHS
HG38FA='/project/holstegelab/Share/pacbio/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa'
CHM13FAGZ='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0.fa.gz'
CHM13FA='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0.fa'
MODELDIR="/project/holstegelab/Share/oscar/software/CpG-main/models/pileup_calling_model"
H38CCS='/project/holstegelab/Share/pacbio/resources/h38_ccs.mmi'
CHM13CCS='/project/holstegelab/Share/asalazar/data/chm13/assembly/v2_0/chm13v2.0_hifi.mmi'
DCACHE_CONFIG='/project/holstegelab/Data/dcache.conf'

CACHE=None

### MAIN -- run before the rules -- identify input files and do md5sum check of the files
inp_bam_prefix = config["IN_FILES"].split(',')
out_bam_prefix = config["OUT_FILE"]
sample_type = out_bam_prefix.split('/')[-2]
dcache_folder = ""
files_to_copy = ""
files_to_copy_names = []
if sample_type in ['ad', 'centenarian']:
    dcache_folder = 'ad_centenarians'
    for s in inp_bam_prefix:
        tmp_files = os.popen('ls ' + s + '.ccs.*').read().replace('\n', ',')
        print(tmp_files)
        fnames = [x.split('/')[-1] for x in tmp_files.split(',') if x.split('/')[-1] != ""]
        files_to_copy = files_to_copy + tmp_files
        files_to_copy_names = files_to_copy_names + fnames

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
        # flye assembly and alignment
        expand("{out_prefix}.hifi.flye/assembly_hg38.bam", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/assembly_chm13.bam", out_prefix = out_bam_prefix),

        # 6. busco analyses
        expand("{out_prefix}_primary_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_alternate_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap1_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}_hap2_BUSCO/logs/busco.log", out_prefix = out_bam_prefix),
        expand("{out_prefix}.hifi.flye/{sample_name}_BUSCO/logs/busco.log", out_prefix = out_bam_prefix, sample_name = out_bam_prefix.split('/')[-1]),

        # 7. SV calling
        # pbsv
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.svcall.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.svcall.vcf'),
        # sniffles
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.hg38.sniffles.snf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.vcf'),
        expand("{inp_prefix}", inp_prefix = out_bam_prefix + '.merged.hifi.chm13.sniffles.snf'),

        # Clean up data from single smrt cells
        ['/project/holstegelab/Software/snakemake_pipeline/trashbin/' + x for x in files_to_copy_names]


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
        {PYTHON39} {PILEUP} -b {input[0]} -f {CHM13FAGZ} -o {params.pfx} -t 50 -p model -d {MODELDIR} --skip_bw_conversion
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
        {SH} {DEEPVARIANT} {input[0]} {output[0]} {output[1]} {CHM13FAGZ}
        """

# Rule to convert BAM file to FASTA file
rule bam2fasta:
    input:
        expand("{out_prefix}.merged.hifi.hg38.bam", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.merged.hifi.fasta", out_prefix = out_bam_prefix))
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

# Rule to align flye assembly -- hg38
rule align_flye_hg38:
    input:
        expand("{out_prefix}.hifi.flye/assembly.fasta", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.flye/assembly_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        """

# Rule to align flye assembly -- chm13
rule align_flye_chm13:
    input:
        expand("{out_prefix}.hifi.flye/assembly.fasta", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.flye/assembly_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        """

# Rule to do BUSCO analysis on flye assembly
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

# Alignment to GRCh38 of hifiasm haplotypes
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

# Alignment to GRCh38 of hifiasm primary
rule align_hifiasm_hg38_primary:
    input:
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.p_ctg_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        """

# Alignment to GRCh38 of hifiasm alternate
rule align_hifiasm_hg38_alternate:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.a_ctg_hg38.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {HG38FA} {input[0]} > {output[0]}
        """

# Alignment to CHM13 of hifiasm haplotypes
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

# Alignment to CHM13 of hifiasm primary
rule align_hifiasm_chm13_primary:
    input:
        expand("{out_prefix}.hifi.hifiasm.p_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.p_ctg_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        """

# Alignment to CHM13 of hifiasm alternate
rule align_hifiasm_chm13_alternate:
    input:
        expand("{out_prefix}.hifi.hifiasm.a_ctg.fa", out_prefix = out_bam_prefix)
    output:
        temp(expand("{out_prefix}.hifi.hifiasm.a_ctg_chm13.sam", out_prefix = out_bam_prefix))
    shell: """
        {MINIMAP2} -ax asm10 -t 20 {CHM13FA} {input[0]} > {output[0]}
        """

# Convert SAM to BAM
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

# SV callling with pbsv -- step 1: discover SV signatures
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

# SV calling with pbsv -- step 2: calling SVs
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

# SV calling with sniffles
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

# Copy single-smrt cells data to dcache when the merging is done and remove them from storage
rule rclone_copy_check_move:
    input:
        [x for x in files_to_copy.split(',') if x != '']
    output:
        ['/project/holstegelab/Software/snakemake_pipeline/trashbin/' + x for x in files_to_copy_names]
    run:
        for f in input:
            print(f)
            # rclone copy
            os.system('rclone copy -vv --progress --multi-thread-streams 1 --config ' + DCACHE_CONFIG + ' ' + f + ' dcache:tape/data_processed/ccs/' + dcache_folder + '/')
            # rclone check
            tmp_check = 'tmp_check_' + str(random.randint(1, 1000000)) + '.txt'
            os.system('rclone check ' + f + ' ~/dcache/tape/data_processed/ccs/' + dcache_folder + '/ --combined ' + tmp_check + ' --size-only --config ' + DCACHE_CONFIG)
            check_res = os.popen("awk '{print $1}' " + tmp_check).read().replace('\n', '')
            if check_res == '=':
                os.system('rm ' + tmp_check)
                os.system('mv ' + f + ' /project/holstegelab/Software/snakemake_pipeline/trashbin/')
