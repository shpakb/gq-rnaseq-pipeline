# include: "raw_data_seq.smk"
# include: "raw_data_chip.smk"
# include: "gds_check.smk"
# include: "gq_processing.smk"
# include: "pca_processing.smk"

configfile: "config.yaml"

import pandas as pd

rule all:
    input:
        expand("out/{organism}/seq/data/filtering/postquant/gsm_stats.tsv",
            organism=["hs", "mm", "rn"])

rule series_matrices_download:
    '''
    Takes search output .txt, parses out all links and downloads them. (wget -nc) Existing up to date files out dir
    doesn't reload.
    '''
    input:
        "input/{organism}/{platform}/gds_search_result.txt"
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        directory("out/{organism}/{platform}/series_matrices"),
    message: "Downloading series matrices for RNA-seq..."
    log: "logs/{organism}/{platform}/series_matrices_seq_download.log"
    shell:
        "scripts/bash/download_sm.sh {input} {output}"
        " > {log} 2>&1"

rule sm_seq_metadata:
    input:
        rules.series_matrices_download.output
    output:
        gsm_df="out/{organism}/{platform}/data/metadata/gsm.tsv",
        gse_df="out/{organism}/{platform}/data/metadata/gse.tsv"
    message: "Aggregating metadata from series matrices..."
    log: "logs/{organism}/{platform}/sm_seq_metadata.log"
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
        gse_df=rules.sm_seq_metadata.output.gse_df,
        gsm_df=rules.sm_seq_metadata.output.gsm_df,
        gpl_df="input/gpl.tsv",
        priority_gse_list="input/{organism}/{platform}/priority_gse.list"
    params:
        min_spots=lambda wildcards: config["min_spots"][wildcards.organism],
        max_spots=lambda wildcards: config["max_spots"][wildcards.organism],
        organism=lambda wildcards: config["organism"][wildcards.organism]
    output:
        gsm_gse_df="out/{organism}/{platform}/data/filtering/prequant/gsm_gse.tsv",
        gsm_filtering_df="out/{organism}/{platform}/data/filtering/prequant/gsm_filtering.tsv",
        passing_gsm_list="out/{organism}/{platform}/data/filtering/prequant/passing_gsm.list",
        srr_gsm_df="out/{organism}/{platform}/data/filtering/prequant/srr_gsm.tsv"
    message: "Pre-quantification filtering..."
    log: "logs/{organism}/{platform}/prequant_filter.log"
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
    srr_df_file = checkpoints.prequant_filter.get(**wildcards).output.srr_gsm_df
    srr_df = pd.read_csv(srr_df_file, sep="\t")
    srr_list = srr_df[srr_df['GSM']==wildcards.gsm]["SRR"].tolist()
    srr_files = \
        expand("out/{organism}/{platform}/kallisto/{srr}/abundance.tsv",
        srr=srr_list,
        organism=wildcards.organism,
        platform=wildcards.platform)
    return srr_files

rule srr_to_gsm:
    resources:
        mem_ram=1
    input:
        srr_files=get_srr_files,
        srr_gsm_df=lambda wildcards: checkpoints.prequant_filter.get(**wildcards).output.srr_gsm_df,
        transcript_gene="input/{organism}/{platform}/transcript_gene.tsv"
    output:
        "out/{organism}/{platform}/gsms/{gsm}.tsv"
    log: "logs/{organism}/{platform}/srr_to_gsm/{gsm}.log"
    message: "Aggregating {wildcards.gsm}"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/srr_to_gsm.R {output} {input.transcript_gene} {input.srr_gsm_df}"
        " > {log} 2>&1"

def get_prequant_filtered_gsm_files(wildcards):
    filtered_gsm_list = checkpoints.prequant_filter.get(**wildcards).output.passing_gsm_list
    gsm_list = [line.rstrip('\n') for line in open(filtered_gsm_list)]
    gsm_files=\
        expand("out/{organism}/{platform}/gsms/{gsm}.tsv",
            gsm=gsm_list,
            organism=wildcards.organism,
            platform=wildcards.platform)
    return gsm_files

checkpoint postquant_filter:
    input:
        gsm_files=get_prequant_filtered_gsm_files,
        gsm_gse_df=lambda wildcards: checkpoints.prequant_filter.get(**wildcards).output.gsm_gse_df
    output:
        gsm_stats_df="out/{organism}/{platform}/data/filtering/postquant/gsm_stats.tsv",
        gsm_gse_df="out/{organism}/{platform}/data/filtering/postquant/gsm_gse.tsv",
        passing_gse_list="out/{organism}/{platform}/data/filtering/postquant/passing_gse.list"
    message: "Post quantification filtering."
    params:
        min_exp_genes=lambda wildcards: config["min_exp_genes"][wildcards.organism]
    log: "logs/{organism}/{platform}/postquant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/postquant_filter.R {config[quant_min_gsm]} {params.min_exp_genes} {input.gsm_gse_df}"
        " {output.gsm_stats_df} {output.gsm_gse_df} {output.passing_gse_list} {input.gsm_files}"
        " > {log} 2>&1"
