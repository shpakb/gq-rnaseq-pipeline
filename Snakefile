# TODO 1: add GSE filtration to the pipeline and exclude all the superserieses
# TODO 2: ...
configfile: "configs/config_rn.yaml"

import pandas as pd
from os import path

#GSES = [line.rstrip('\n') for line in open(config["gse_list"])]
GSES = ["GSE61179"]

SRR_GSM_DF = pd.read_csv(config["srr_to_gsm"], sep="\t")
GSM_GSE_DF = pd.read_csv(config["gsm_to_gse"], sep="\t")

rule all:
    input:
        "out/stats/gsm_stats.tsv",
        "out/stats/gse_gsm_filtered.tsv",
        "out/flag"

rule sra_download:
    resources:
        download_res=25,
        mem_ram=1
    output: temp("out/sra/{srr}.sra")
    log:    "out/sra/{srr}.log"
    message: "Downloading {wildcards.srr}"
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/download_sra.sh {wildcards.srr} {output} > {log} 2>&1"

# writes complete flag file to output folder
# consumes writing resources
rule sra_fastqdump:
    resources:
        writing_res=1
    input:
        "out/sra/{srr}.sra"
    output:
        fastq_dir=temp(directory("out/fastq/{srr}")),
        complete_flag="out/fastq/{srr}_complete"
    log:    "out/fastq/{srr}.log"
    message: "fastq-dump {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shell:
        "fastq-dump --outdir {output.fastq_dir} --split-3 {input} >{log} 2>&1 &&"
        "touch {output.complete_flag}"

rule fastq_kallisto:
    resources:
        mem_ram=4
    input:
        rules.sra_fastqdump.output.complete_flag,
        fastq_dir="out/fastq/{srr}"
    output:
        h5=protected("out/kallisto/{srr}/abundance.h5"),
        tsv=protected("out/kallisto/{srr}/abundance.tsv"),
        json=protected("out/kallisto/{srr}/run_info.json")
    log: "out/kallisto/{srr}/{srr}.log"
    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/quantify.sh {wildcards.srr} {input.fastq_dir} "
        " {config[refseq]} out/kallisto/{wildcards.srr} >{log} 2>&1"

rule srr_to_gsm:
    resources:
        mem_ram=1
    input:
        lambda wildcards: expand(
            "out/kallisto/{srr}/abundance.tsv",
            srr=SRR_GSM_DF[SRR_GSM_DF.gsm == wildcards.gsm]["srr"].tolist()
        )
    output: "out/gsms/{gsm}.tsv"
    log: "out/gsms/{gsm}.log"
    message: "Aggregating GSM {wildcards.gsm}"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/srr_to_gsm.R {wildcards.gsm}"
        " {config[probes_to_genes]} {config[srr_to_gsm]}"
        " out/kallisto out/gsms"

checkpoint get_gsm_stats:
    resources:
        mem_ram=2
    input:
        gsm_files=lambda wildcards: expand(
            "out/gsms/{gsm}.tsv",
            gsm=GSM_GSE_DF[GSM_GSE_DF['gse'].isin(GSES)]["gsm"].tolist()
        ),
    output:
        gsm_stats="out/stats/gsm_stats.tsv",
        gse_filtered_df="out/stats/gse_gsm_filtered.tsv"
    shell:
        "Rscript scripts/get_gsm_stats.R {output.gsm_stats} {input.gsm_files} &&"
        "   Rscript scripts/filtered_gse_list.R {output.gsm_stats}"
        "       {output.gse_filtered_df} {config[gsm_to_gse]} {config[min_gse_gq]} {config[min_exp_genes]}"


# def get_gsm_list(wildcards):
#     gse_filtered_path = checkpoints.get_gsm_stats.get(**wildcards).output.gse_filtered_df
#     gse_gsm_df = pd.read_csv(gse_filtered_path, sep="\t")
#     return gse_gsm_df[gse_gsm_df.gse == wildcards.gse]["gsm"].tolist()


rule gsm_to_gse:
    input:
        gse_filtered_df="out/stats/gse_gsm_filtered.tsv"
    output:
        gse="out/gses/{gse}.tsv"
    log: "out/gses/{gse}.log"
    message: "Aggregating GSE"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/gsm_to_gse.R {output} out/gsms"
        " {config[ensamble_genesymbol_entrez]} "
        " {input.gse_filtered_df}"

def get_filtered_gses(wildcards):
    gse_filtered_path = checkpoints.get_gsm_stats.get(**wildcards).output.gse_filtered_df
    gse_gsm_df = pd.read_csv(gse_filtered_path, sep="\t")
    gses = list(set(gse_gsm_df["gse"].tolist()))
    gse_files = ["out/gses/" + gse + ".tsv" for gse in gses]
    return gse_files

# rule for pushing filtered gses through the pipeline
rule push_filtered_gses:
    input:
        get_filtered_gses
    output:
        flag="out/flag"
    shell:
        "touch {output.flag}"

# rule cluster:
#     output:
#         "out/wgcna/{gse}/processing.log"
#     shell:
#         "Rscript scripts/wgcna_seq.R"
