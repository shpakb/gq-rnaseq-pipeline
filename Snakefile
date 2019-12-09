configfile: "config.yaml"

# TODO: PCA query on rna-seq
# TODO: aggregate metadata on kallisto runs in df and choose version of callisto used for most of the samples.


import pandas as pd
from os import listdir
from os.path import isfile, join
import re

rule all:
    input:
        # expand("out/{organism}/chip/sm_qc/qc_df.tsv",
        #     organism=["hs", "mm", "rn"]),
        # expand("out/{organism}/seq/prequant_filter/gsm_gse.tsv",
        #     organism=["hs", "mm", "rn"]),
        # expand("out/{organism}/{platform}/sm_metadata/gsm.tsv",
        #     organism=["hs", "mm", "rn"],
        #     platform=["chip", "seq"]),
        # "out/data/srr_gsm_spots.tsv",
        expand("out/{organism}/chip/exp_qc_df.tsv",
            organism=["rn"])
        # expand("out/{organism}/{platform}/wgcna_stats.tsv",
        #     organism=["hs", "mm", "rn"],
        #     platform=["chip", "seq"])

rule sm_download:
    '''
    Takes search output .txt, parses out all links and downloads them. (wget -nc) Existing up to date files out dir
    doesn't reload. Writes completion flag in the end.
    '''
    resources:
        time=60*24*2,
        download_res=1,
        writing_res=1,
        priority=5
    input:
        "input/{organism}/{platform}/gds_search_result.txt"
    output:
        "flags/{organism}/{platform}/sm_download.flag"
    params:
        sm_download_dir="out/{organism}/{platform}/series_matrices"
    message: "Downloading series matrices for {wildcards.organism} {wildcards.platform}..."
    log: "logs/{organism}/{platform}/sm_download.log"
    shell:
        "scripts/bash/download_sm.sh {input} {params.sm_download_dir} "
        " > {log} 2>&1 &&"
        " touch {output}"

checkpoint extract_sm_metadata:
    input:
        rules.sm_download.output
    output:
        gsm_df="out/{organism}/{platform}/sm_metadata/gsm.tsv",
        gse_df="out/{organism}/{platform}/sm_metadata/gse.tsv"
    params:
        sm_download_dir=rules.sm_download.params.sm_download_dir
    message: "Aggregating metadata from series matrices {wildcards.organism} {wildcards.platform} ..."
    log: "logs/{organism}/{platform}/sm_seq_metadata.log"
    shell:
        "python scripts/python/parse_sm_metadata.py {params.sm_download_dir} {output.gse_df} {output.gsm_df}"
        " > {log} 2>&1"

########################################################################################################################
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
    message: "Pre-quantification filtering for {wildcards.organism}..."
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
    message: "Downloading {wildcards.srr} ({wildcards.organism})..."
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
    message: "fastq-dump {wildcards.srr} ({wildcards.organism})"
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
    message: "Kallisto: {wildcards.srr} ({wildcards.organism})"
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
    message: "Aggregating {wildcards.gsm} ({wildcards.organism})..."
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
    message: "Post quantification filtering ({wildcards.organism}) ..."
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
    message: "Aggregating {wildcards.gse} ({wildcards.organism})..."
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
        expand("out/{organism}/seq/gses/{gse}.tsv",
            gse=gses,
            organism=wildcards.organism)
    return gse_files

rule push_gse:
    input: get_postquant_passing_gse
    output: "out/{organism}/seq/push_gse_flag"
    shell: "touch {output}"

################################################CHIP_ROOT###############################################################
rule extract_exp_mat:
    '''
    Extracts expression matrices from files and writes QC report for downstream filtering.
    qc_df: TAG N_GSM HAS_EXP_MAT GPL LOGAV LOGMAX LINMAX N_GENES
    Writes empty exp table if error produced.
    '''
    input:
        sm="out/{organism}/chip/series_matrices/{tag}_series_matrix.txt.gz",
        gpl_dir="input/{organism}/chip/platform_annotation"
    output:
        exp_table="out/{organism}/chip/exp_mat/{tag}.tsv",
        qc_report="out/{organism}/chip/exp_mat/{tag}_qc.tsv"
    log: "logs/{organism}/chip/extract_exp_mat/{tag}.log"
    message: "Extracting expression table from {wildcards.tag}_series_matrix.txt.gz ({wildcards.organism})..."
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/get_exp_table.R {input.sm} {input.gpl_dir} {output.exp_table} {output.qc_report}"
        " > {log} 2>&1"


def get_filtered_sm_qc_files(wildcards):
    gse_df_file = checkpoints.extract_sm_metadata.get(**wildcards).output.gse_df
    gse_df = pd.read_csv(gse_df_file, sep="\t")
    gse_df = gse_df[(gse_df['NUMBER_GSM']>=config['min_gsm']) | (gse_df['NUMBER_GSM']<=config['max_gsm'])]
    # getting list of supported GPL
    gpl_dir = f"input/{wildcards.organism}/chip/platform_annotation"
    gpl_list = [f for f in listdir(gpl_dir) if isfile(join(gpl_dir, f))]
    gpl_list = [re.search('(.+?).3col.gz', gpl).group(1) for gpl in gpl_list]
    gse_df = gse_df[gse_df["GPL"].isin(gpl_list)]

    gse_df = gse_df[gse_df["SUPER_SERIES_OF"].isna()]

    filtered_sm_tags=gse_df["GSE"].to_list()

    filtered_sm_qc_files = \
        expand('out/{organism}/chip/exp_mat/{tag}_qc.tsv',
        organism=wildcards.organism,
        tag=filtered_sm_tags)

    return filtered_sm_qc_files

rule get_exp_mat_qc_df:
    input:
        get_filtered_sm_qc_files
    output:
        "out/{organism}/{platform}/exp_qc_df.tsv"
    message: "Getting expression matrices QC df ({wildcards.organism})..."
    log: "logs/{organism}/{platform}/get_exp_mat_qc_df.log"
    run:
        qc_df_list=[pd.read_csv(qc_file,sep='\t') for qc_file in input]
        qc_df=pd.concat(qc_df_list)
        qc_df.to_csv(str(output),sep='\t',index=False)
        print("Done.")


###############################################WGCNA####################################################################
def get_exp_mat_for_wgcna(wildcards):
    """
    Some qc logic and rooting for choosing files for downstream WGCNA.
    Writes folder with 3 files. Modules, eigengenes and WGCNA stats. If WGCNA failed, Modules and Iegenegenes files are
    empty and stats contains reason WGCNA failed.
    params:
        logav_min=config["chip_logav_min"],
        logav_max=config["chip_logav_max"],
        linmax_max=config["chip_linmax_max"],
        logmax_max=config["chip_logmax_max"],
        min_genes=config["min_genes"]
    """
    if wildcards.platform=="chip":
        sm_qc_df_file = rules.extract_exp_mat.output.sm_qc_df # WON'T WORK
        sm_qc_df = pd.read_csv(sm_qc_df_file, sep="\t")
        sm_qc_df = sm_qc_df[sm_qc_df['PROCESSED']==True]
        sm_qc_df = sm_qc_df[sm_qc_df['N_GSM']>=config['min_gsm_gq'] | sm_qc_df['N_GSM']<=config['max_gsm_gq']]
        exp_mat_tags = sm_qc_df["TAG"].tolist()
        exp_mat_files = \
            expand("out/{organism}/chip/exp_mat/{tag}.tsv",
                tag=exp_mat_tags,
                organism=wildcards.organism)
    elif wildcards.platform=="seq":
        gsm_gse_df_file = checkpoints.postquant_filter.get(**wildcards).output.gsm_gse_df
        gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
        gse_df =  gse_df.groupby('GSE')['GSE'].transform('count')
        exp_mat_list = gse_df[gse_df["freq"]>=config['min_gsm_gq'] | gse_df['freq']<=config['max_gsm_gq']].tolist()
        exp_mat_files = \
            expand("out/{organism}/chip/gses/{tag}.tsv",
            tag=exp_mat_list,
            organism=wildcards.organism)

    return exp_mat_files

rule wgcna:
    input:
        get_exp_mat_for_wgcna
    output:
        modules="out/{organism}/{platform}/wgcna/{tag}/modules.tsv",
        eigengenes="out/{organism}/{platform}/wgcna/{tag}/eigengenes.tsv",
        stats="out/{organism}/{platform}/wgnca/{tag}/stats.txt"
    params:
        n_genes=config["wgcna_n_genes"],
        cor_type=config["wgcna_cor"],
        network_type=config["wgcna_nt"],
        logscale_f=config["wgcna_logscale"]
    log: "logs/{organism}/{platform}/wgcna/{tag}.log"
    message: "Performing WGCNA on {wildcards.tag} ({wildcards.organism}, {wildcards.platform})..."
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/wgcna.R {input} {output.modules} {output.eigengenes} {output.stats} "
        " {params.n_genes} {params.cor_type} {params.network_type} {params.logscale_f}"
        " > {log} 2>&1"

checkpoint get_wgcna_stats:
    output:
        "out/{organism}/{platform}/wgcna_stats.tsv"
    params:
        "out/{organism}/{platform}/wgcna"
    conda: "envs/r_scripts.yaml"
    log: "logs/{organism}/{platform}/get_wgcna_stats.log"
    message: "Getting WGCNA stats ({wildcards.organism}, {wildcards.platform})..."
    shell:
        "Rscript scripts/R/get_wgcna_stats.R {params} {output}"
        " > {log} 2>&1"

