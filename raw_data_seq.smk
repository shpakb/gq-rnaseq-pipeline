# Issues: can't pass files by input functions from checkpoints in runtime. But still works during graph calculation.
# 1) Cant pass temp(directory("out/fasqdump/")) as reference to from rules variables. Had to make workaround with
# completion flags

# Add list with priority GSE that bypasses both checkpoints
# (except for the filtering on right organism and right type of samples)

# TODO: move all QC logic to input functions(if possible)

import pandas as pd

configfile: "configs/config_rn.yaml"

rule series_matrices_seq_download:
    '''
    Takes search output .txt, parses out all links and downloads them. 
    (wget -nc) Existing up to date files out dir doesn't reload. To update  database just remove complete flag file. 
    '''
    input:
        config["geo_search_result_chip"]
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        sm_dir=directory("out/series_matrices/seq/"),
        complete_flag="out/flags/series_matrices_seq_download"
    message: "Downloading series matrices for RNA-seq..."
    log: "out/logs/series_matrices_seq_download.log"
    shell:
        "scripts/bash/download_sm.sh {input} {output.sm_dir}"
        " > {log} 2>&1 &&"
        " touch {output.complete_flag}"

rule sm_seq_metadata:
    input:
        rules.series_matrices_seq_download.output.sm_dir
    output:
        gsm_df="out/data/metadata/seq/gsm.tsv",
        gse_df="out/data/metadata/seq/gse.tsv"
    message: "Aggregating metadata from series matrices..."
    log: "out/logs/sm_seq_metadata.log"
    shell:
        "python scripts/python/parse_sm_metadata.py {input} {output.gse_df} {output.gsm_df}"
        " > {log} 2>&1"

rule sra_accession_df_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        temp("out/data/sra_accession_raw.tab")
    message: "Downloading SRA accession table..."
    log: "out/logs/sra_accession_df_download.log"
    shell:
        "wget -O {output} {config[sra_accession_df]}"
        " > {log} 2>&1"

rule get_srr_df:
    resources:
        writing_res=1
    input:
        rules.sra_accession_df_download.output
    output:
        srr_df="out/data/srr_gsm_spots.tsv"
    message: "Cleaning SRR df..."
    log: "out/logs/get_srr_df.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "scripts/bash/clean_sra_accession_df.sh {input} {output.srr_df}"
        " > {log} 2>&1 &&"
        " Rscript scripts/R/clean_sra_accession_df.R {output.srr_df}"
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
        srr_df=rules.get_srr_df.output.srr_df,
        gse_df=rules.sm_seq_metadata.output.gse_df,
        gsm_df=rules.sm_seq_metadata.output.gsm_df,
        gpl_df=config["gpl_df"],
        priority_gse_list=config["priority_gse_list"]
    output:
        gsm_gse_df="out/data/filtering/prequant/gsm_gse.tsv",
        gsm_filtering_df="out/data/filtering/prequant/gsm_filtering.tsv",
        passing_gsm_list="out/data/filtering/prequant/passing_gsm.list",
        srr_gsm_df="out/data/filtering/prequant/srr_gsm.tsv"
    message: "Pre-quantification filtering..."
    log: "out/logs/prequant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/prequant_filter.R {input.gse_df} {input.gsm_df} {input.gpl_df} {input.srr_df}"
        " {config[organism]} {config[min_spots_gsm]} {config[max_spots_gsm]} {config[quant_min_gsm]} "
        " {config[quant_max_gsm]} {output.gsm_filtering_df} {output.passing_gsm_list} {output.srr_gsm_df}"
        " {output.gsm_gse_df} {input.priority_gse_list} {config[priority_only_f]}"
        " > {log} 2>&1"

rule sra_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output: temp("out/sra/{srr}.sra")
    log:    "out/logs/sra_download/{srr}.log"
    message: "Downloading {wildcards.srr}"
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/bash/download_sra.sh {wildcards.srr} {output}"
        " > {log} 2>&1"

# TODO: migrate to fasterq-dump
rule sra_fastqdump:
    input:
        "out/sra/{srr}.sra"
    output:
        fastq_dir=temp(directory("out/fastq/{srr}")),
        complete_flag=temp("out/fastq/{srr}_complete")
    log:    "out/logs/sra_fastqdump/{srr}.log"
    message: "fastq-dump {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shell:
        "fastq-dump --outdir {output.fastq_dir} --split-3 {input}"
        " > {log} 2>&1 &&"
        " touch {output.complete_flag}"

rule fastq_kallisto:
    resources:
        mem_ram=8
    input:
        rules.sra_fastqdump.output.complete_flag,
        fastq_dir="out/fastq/{srr}"
    output:
        h5=protected("out/kallisto/{srr}/abundance.h5"),
        tsv=protected("out/kallisto/{srr}/abundance.tsv"),
        json=protected("out/kallisto/{srr}/run_info.json")
    log: "out/logs/fastq_kallisto/{srr}.log"
    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/bash/quantify.sh {wildcards.srr} {input.fastq_dir} {config[refseq]} out/kallisto/{wildcards.srr}"
        " > {log} 2>&1"

# For some reason input function returns proper output during graph calculation but
# returns first output variable from checkpoint during rule execution
# This does guaranties that all necessary files for the task are there but doesn't allow to pass them as an input to
# script.


def get_srr_files(wildcards):
    srr_df_file = checkpoints.prequant_filter.get(**wildcards).output.srr_gsm_df
    srr_df = pd.read_csv(srr_df_file, sep="\t")
    srr_list = srr_df[srr_df['GSM']==wildcards.gsm]["SRR"].tolist()
    srr_files = expand("out/kallisto/{srr}/abundance.tsv", srr=srr_list)
    return srr_files

rule srr_to_gsm:
    resources:
        mem_ram=1
    input:
        srr_files=get_srr_files,
        srr_df="out/data/filtering/prequant/srr_gsm.tsv"
    output:
        gsm_file="out/gsms/{gsm}.tsv"
    log: "out/logs/srr_to_gsm/{gsm}.log"
    message: "Aggregating {wildcards.gsm}"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/srr_to_gsm.R {output.gsm_file} {config[probes_to_genes]} {input.srr_df}"
        " > {log} 2>&1"

def get_prequant_filtered_gsm(wildcards):
    filtered_gsm_list = checkpoints.prequant_filter.get(**wildcards).output.passing_gsm_list
    gsms = [line.rstrip('\n') for line in open(filtered_gsm_list)]
    gsm_files = expand("out/gsms/{gsm}.tsv", gsm=gsms)
    return gsm_files

def get_gsm_gse_df(wildcards):
    return checkpoints.prequant_filter.get(**wildcards).output.gsm_gse_df

checkpoint postquant_filter:
    input:
        gsm_files=get_prequant_filtered_gsm,
        gsm_gse_df=get_gsm_gse_df
    output:
        gsm_stats_df="out/data/filtering/postquant/gsm_stats.tsv",
        gsm_gse_df="out/data/filtering/postquant/gsm_gse.tsv",
        passing_gse_list="out/data/filtering/postquant/passing_gse.list"
    message: "Post quantification filtering."
    log: "out/logs/postquant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/postquant_filter.R {config[quant_min_gsm]} {config[min_exp_genes]} {input.gsm_gse_df}"
        " {output.gsm_stats_df} {output.gsm_gse_df} {output.passing_gse_list} {input.gsm_files}"
        " > {log} 2>&1"

def get_postquant_gsms_for_gse(wildcards):
    gsm_gse_df_file = checkpoints.postquant_filter.get(**wildcards).output.gsm_gse_df
    gsm_gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
    gsms = gsm_gse_df[gsm_gse_df['GSE']==wildcards.gse]["GSM"].tolist()
    gsm_files = expand("out/gsms/{gsm}.tsv", gsm=gsms)
    return gsm_files

rule gsm_to_gse:
    input:
        gsm_files=get_postquant_gsms_for_gse,
        gse_df="out/data/filtering/postquant/gsm_gse.tsv"
    output:
        gse="out/gses/{gse}.tsv"
    log: "out/logs/gsm_to_gse/{gse}.log"
    message: "Aggregating {wildcards.gse}"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/gsm_to_gse.R {output.gse} out/gsms {config[ensamble_genesymbol_entrez]} {input.gse_df}"
        " > {log} 2>&1"

def get_postquant_passing_gse(wildcards):
    filtered_gse_list = checkpoints.postquant_filter.get(**wildcards).output.passing_gse_list
    gses = [line.rstrip('\n') for line in open(filtered_gse_list)]
    gse_files = expand("out/gses/{gsm}.tsv", gsm=gses)
    return gse_files

rule push_filtered_gses:
    input:
        get_postquant_passing_gse
    output:
        flag="out/flag"
    shell:
        "touch {output.flag}"
