configfile: "config.yaml"

import pandas as pd

rule all:
    input:
        expand("out/{organism}/seq/prequant_filter/gsm_gse.tsv",
            organism=["hs", "mm", "rn"]),
        expand("out/{organism}/{platform}/sm_metadata/gsm.tsv",
            organism=["hs", "mm", "rn"],
            platform=["chip", "seq"]),
        "out/data/srr_gsm_spots.tsv"
        # expand("out/{organism}/{platform}/sm_metadata/gsm.tsv",
        #     organism=["rn"],
        #     platform=["chip", "seq"])


rule sm_download:
    '''
    Takes search output .txt, parses out all links and downloads them. (wget -nc) Existing up to date files out dir
    doesn't reload. Writes completion flag in the end.
    '''
    input:
        "input/{organism}/{platform}/gds_search_result.txt"
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        "flags/{organism}/{platform}/sm_download.flag"
    params:
        sm_download_dir="out/{organism}/{platform}/series_matrices"
    message: "Downloading series matrices for RNA-seq..."
    log: "logs/{organism}/{platform}/sm_download.log"
    shell:
        "scripts/bash/download_sm.sh {input} {params.sm_download_dir} "
        " > {log} 2>&1 &&"
        " touch {output}"

rule sm_metadata:
    input:
        rules.sm_download.output
    output:
        gsm_df="out/{organism}/{platform}/sm_metadata/gsm.tsv",
        gse_df="out/{organism}/{platform}/sm_metadata/gse.tsv"
    params:
        sm_download_dir=rules.sm_download.params.sm_download_dir
    message: "Aggregating metadata from series matrices..."
    log: "logs/{organism}/{platform}/sm_seq_metadata.log"
    shell:
        "python scripts/python/parse_sm_metadata.py {params.sm_download_dir} {output.gse_df} {output.gsm_df}"
        " > {log} 2>&1"

rule sra_accession_df_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        temp("out/data/sra_accession_raw.tab")
    message: "Downloading SRA accession table..."
    log: "logs/sra_accession_df_download.log"
    shell:
        "wget -O {output} {config[sra_accession_df]}"
        " > {log} 2>&1"

rule get_srr_df:
    resources:
        writing_res=1
    input:
        rules.sra_accession_df_download.output
    output:
        "out/data/srr_gsm_spots.tsv"
    message: "Cleaning SRR df..."
    log: "logs/get_srr_df.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "scripts/bash/clean_sra_accession_df.sh {input} {output}"
        " > {log} 2>&1 &&"
        " Rscript scripts/R/clean_sra_accession_df.R {output}"
        " >> {log} 2>&1"

checkpoint prequant_filter:
    '''
    Outputs QC tables and lists of passing GSE and GSM. Takes one extra argument: priority_gse_list with GSEs that get
    to downstream quantification regardless of QC. Keep in mind that GSM for priority GSE gets filtered by sample
    metadata(GSM can be microarray, wrong organism etc.) 
    Temporary added priority_only flag and script logic for test runs and to quantify things while pipeline under 
    development.
    '''
    input:
        srr_df=rules.get_srr_df.output,
        gse_df="out/{organism}/seq/sm_metadata/gse.tsv",
        gsm_df="out/{organism}/seq/sm_metadata/gsm.tsv",
        gpl_df="input/gpl.tsv",
        priority_gse_list="input/{organism}/seq/priority_gse.list"
    params:
        min_spots=lambda wildcards: config["min_spots"][wildcards.organism],
        max_spots=lambda wildcards: config["max_spots"][wildcards.organism],
        organism=lambda wildcards: config["organism"][wildcards.organism]
    output:
        gsm_gse_df="out/{organism}/seq/prequant_filter/gsm_gse.tsv",
        gsm_filtering_df="out/{organism}/seq/prequant_filter/gsm_filtering.tsv",
        passing_gsm_list="out/{organism}/seq/prequant_filter/passing_gsm.list",
        srr_gsm_df="out/{organism}/seq/prequant_filter/srr_gsm.tsv"
    message: "Pre-quantification filtering..."
    log: "logs/{organism}/seq/prequant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/prequant_filter.R {input.gse_df} {input.gsm_df} {input.gpl_df} {input.srr_df}"
        " {params.organism} {params.min_spots} {params.max_spots} {config[quant_min_gsm]}"
        " {config[quant_max_gsm]} {output.gsm_filtering_df} {output.passing_gsm_list} {output.srr_gsm_df}"
        " {output.gsm_gse_df} {input.priority_gse_list} {config[priority_only_f]}"
        " > {log} 2>&1"

rule sra_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    priority: 3
    output: temp("out/{organism}/seq/sra/{srr}.sra")
    log:    "logs/{organism}/seq/sra_download/{srr}.log"
    message: "Downloading {wildcards.srr}"
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/bash/download_sra.sh {wildcards.srr} {output}"
        " > {log} 2>&1"

rule sra_fastqdump:
    resources:
        writing_res=1
    input:
        rules.sra_download.output
    output:
        fastq_dir=temp(directory("out/{organism}/seq/fastq/{srr}")),
    log:    "logs/{organism}/seq/sra_fastqdump/{srr}.log"
    message: "fastq-dump {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "fastq-dump --outdir {output.fastq_dir} --split-3 {input}"
        " > {log} 2>&1"

rule fastq_kallisto:
    resources:
        mem_ram=lambda wildcards: config["quant_mem_ram"][wildcards.organism]
    priority: 2
    input:
        fastq_dir=rules.sra_fastqdump.output,
        # fastq_dir="out/{organism}/seq/fastq/{srr}",
        refseq="input/{organism}/seq/refseq_kallisto"
    output:
        h5=protected("out/{organism}/seq/kallisto/{srr}/abundance.h5"),
        tsv=protected("out/{organism}/seq/kallisto/{srr}/abundance.tsv"),
        json=protected("out/{organism}/seq/kallisto/{srr}/run_info.json")
    log: "logs/{organism}/seq/fastq_kallisto/{srr}.log"
    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/bash/quantify.sh {wildcards.srr} {input.fastq_dir} {input.refseq}"
        " out/{wildcards.organism}/seq/kallisto/{wildcards.srr}" # output dir 
        " > {log} 2>&1"

def get_srr_files(wildcards):
    srr_df_file = f"out/{wildcards.organism}/seq/prequant_filter/srr_gsm.tsv"
    srr_df = pd.read_csv(srr_df_file, sep="\t")
    srr_list = srr_df[srr_df['GSM']==wildcards.gsm]["SRR"].tolist()
    srr_files = \
        expand("out/{organism}/seq/kallisto/{srr}/abundance.tsv",
        srr=srr_list,
        organism=wildcards.organism)
    return srr_files

rule srr_to_gsm:
    resources:
        mem_ram=1
    input:
        srr_files=get_srr_files,
        srr_gsm_df="out/{organism}/seq/prequant_filter/srr_gsm.tsv",
        transcript_gene="input/{organism}/seq/transcript_gene.tsv"
    output:
        "out/{organism}/seq/gsms/{gsm}.tsv"
    log: "logs/{organism}/seq/srr_to_gsm/{gsm}.log"
    message: "Aggregating {wildcards.gsm}"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/srr_to_gsm.R {output} {input.transcript_gene} {input.srr_gsm_df}"
        " > {log} 2>&1"

def prequant_filtered_gsm_files(wildcards):
    filtered_gsm_list = checkpoints.prequant_filter.get(**wildcards).output.passing_gsm_list
    gsm_list = [line.rstrip('\n') for line in open(filtered_gsm_list)]
    gsm_files=\
        expand("out/{organism}/seq/gsms/{gsm}.tsv",
            gsm=gsm_list,
            organism=wildcards.organism)
    return gsm_files

checkpoint postquant_filter:
    input:
        gsm_files=prequant_filtered_gsm_files,
        gsm_gse_df=lambda wildcards: checkpoints.prequant_filter.get(**wildcards).output.gsm_gse_df
    output:
        gsm_stats_df="out/{organism}/seq/postquant_filter/gsm_stats.tsv",
        gsm_gse_df="out/{organism}/seq/postquant_filter/gsm_gse.tsv",
        passing_gse_list="out/{organism}/seq/postquant_filter/passing_gse.list"
    message: "Post quantification filtering."
    params:
        min_exp_genes=lambda wildcards: config["min_exp_genes"][wildcards.organism]
    log: "logs/{organism}/seq/postquant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/postquant_filter.R {config[quant_min_gsm]} {params.min_exp_genes} {input.gsm_gse_df}"
        " {output.gsm_stats_df} {output.gsm_gse_df} {output.passing_gse_list} {input.gsm_files}"
        " > {log} 2>&1"

def postquant_gsms_for_gse(wildcards):
    gsm_gse_df_file = f"out/{wildcards.organism}/seq/postquant_filter/gsm_gse.tsv"
    gsm_gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
    gsms = gsm_gse_df[gsm_gse_df['GSE']==wildcards.gse]["GSM"].tolist()
    gsm_files = \
        expand("out/{organism}/seq/gsms/{gsm}.tsv",
            gsm=gsms,
            organism=wildcards.organism)
    return gsm_files

rule gsm_to_gse:
    input:
        gsm_files=postquant_gsms_for_gse,
        gse_df="out/{organism}/seq/data/filtering/postquant/gsm_gse.tsv",
        gene_mapping="input/{organism}/ensamble_genesymbol_entrez.tsv}"
    output:
        gse="out/{organism}/seq/gses/{gse}.tsv"
    log: "logs/{organism}/gsm_to_gse/{gse}.log"
    message: "Aggregating {wildcards.gse}"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/gsm_to_gse.R {output.gse} out/{wildcards.organism}/seq/gsms {input.gene_mapping} "
        " {input.gse_df}"
        " > {log} 2>&1"

def get_postquant_passing_gse(wildcards):
    filtered_gse_list = checkpoints.postquant_filter.get(**wildcards).output.passing_gse_list
    gses = [line.rstrip('\n') for line in open(filtered_gse_list)]
    gse_files = \
        expand("out/{organism}/seq/gses/{gsm}.tsv",
            gsm=gses,
            organism=wildcards.organism)
    return gse_files

rule push_gse:
    input: get_postquant_passing_gse
    output: "out/{organism}/seq/push_gse_flag"
    shell: "touch {output}"

################################################CHIP_ROOT###############################################################
checkpoint extract_exp_mat:
    '''
    Takes SM folder and gse.tsv. 
    Filters out SM by number of GSM, GPL. If passing first two, enters the SM dir and opens the corresponded SM file. 
    Checks: 
        1) If exp table is present
        2) If it can be read properly 
        3) If passes logAverage linAverage QC (might be something else) 
    Outputs sm qc_df where columns contain values for QC and the last column PASSING
    '''
    input:
        gse_df=rules.sm_metadata.output.gse_df,
        gpl_dir="input/{organism}/chip/3col"
    params:
        sm_download_dir=rules.sm_download.params.sm_download_dir
    output:
        sm_qc_df="out/{organism}/chip/prefilter_chip_sm/qc_df.tsv",
        qc_values_file="out/{organism}/chip/prefilter_chip_sm/qc_params.txt",
        exp_df_dir=directory("out/{organism}/chip/exp_mat")
    log: "logs/{organism}/chip/extract_exp_mat.log"

    shell:
        "Rscript chip_qc.R {input.gse_df} {params.sm_download_dir} {input.gpl_dir} {output.sm_qc_df}"
        " {output.qc_values_file} config[min_gsm_chip] config[max_gsm_chip]"
        " > {log} 2>&1"

def prefiltered_chip_sm(wildcards):
    sm_qc_df_file=checkpoints.extract_exp_mat.get(**wildcards).output.sm_qc_df
    sm_qc_df = pd.read_csv(sm_qc_df_file, sep="\t")
    exp_mat_tags = sm_qc_df[sm_qc_df['PASSED']==True]["TAG"].tolist()
    exp_mat_files = \
        expand("out/{organism}/chip/exp_mat/{tag}.tsv",
            tag=exp_mat_tags,
            organism=wildcards.organism)
    return exp_mat_files

rule get_exp_matrix:
    '''
    Trues to extract expression matrix from SM. Performs QC. If fail to extract matrix writes failure reason right to
    the expression table file. 
    '''
    input: "input"
    output: "output"

checkpoint finalize_exp_tables:
    '''
    Gets 
    '''
