# TODO 1: add GSE filtration to the pipeline and exclude all the superserieses
# TODO 2: ...
configfile: "configs/config_rn.yaml"

import pandas as pd
# import os

#GSES = [line.rstrip('\n') for line in open(config["gse_list"])]
GSES = ["GSE61179"]

SRR_GSM_DF = pd.read_csv(config["srr_to_gsm"], sep="\t")
GSM_GSE_DF = pd.read_csv(config["gsm_to_gse"], sep="\t")

rule all:
    input:
        "out/stats/gsm_stats.tsv",
        "out/stats/gse_gsm_filtered.tsv",
        "out/flags/gses_to_merge_collected"

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
        "scripts/download_sra.sh {wildcards.srr} {output} >{log} 2>&1"

# TODO: output empty folder in case of failure
rule sra_kallisto_quant:
    resources:
        sra_kallisto_quant=1,
        mem_ram=4
    input:
        sra="out/sra/{srr}.sra",
    output:
        h5=protected("out/kallisto/{srr}/abundance.h5"),
        tsv=protected("out/kallisto/{srr}/abundance.tsv"),
        json=protected("out/kallisto/{srr}/run_info.json")
    log: "out/kallisto/{srr}/{srr}.log"
    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/quantify.sh {wildcards.srr} {input.sra} "
        " {config[refseq]} out/kallisto/{wildcards.srr} >{log} 2>&1"

# TODO: sepparate on two fastq-split and kallisto quant
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

# checkpoint rule that creates df:
# GSM_ID N_GENES_EXP. If GSM file empty N_GENES_EXP is NA
# And outputs list of GSE to aggregate
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

def collect_gses_to_merge_input_fun(wildcards):
    gse_filtered_path = checkpoints.get_gsm_stats.get(**wildcards).output.gse_filtered_df
    gse_gsm_df = pd.read_csv(gse_filtered_path, sep="\t")
    return [f"out/gses/{gse}.tsv" for gse in gse_gsm_df['gse'].tolist()]

rule collect_gses_to_merge:
    input: collect_gses_to_merge_input_fun
    output: touch("out/flags/gses_to_merge_collected")

def get_gsm(gse_id, gse_gsm_tsv):
    df = pd.read_csv(gse_gsm_tsv, sep="\t")
    res = df[df['gse'] == gse_id]["gsm"].tolist()
    return res


rule gsm_to_gse:
    input:
        gsms=lambda wildcards: expand(
            "out/gsms/{gsm}.tsv",
            gsm=get_gsm({wildcards.gse}, checkpoints.get_gsm_stats.get(**wildcards).output.gse_filtered_df)
        )
    output: "out/gses/{gse}.tsv"
    log: "out/gses/{gse}.log"
    message: "Aggregating GSE"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/gsm_to_gse.R {wildcards.gse}"
        " {config[ensamble_genesymbol_entrez]} "
        " out/stats/gse_gsm_filtered.tsv {input.gsms}"

# checkpoint gse_stats:
#     input:
#         expand(
#             "out/gses/{gse}.tsv",
#             gse=GSES
#         )
#     output:
#         "out/gse_stats.tsv"
#     run:
#         gse_passed_gq = []
#         for gse_f in input:
#             with open(gse_f) as f:
#                 reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
#                 first_row = next(reader)
#                 num_cols = len(first_row)
#                 if os.path.getsize(gse_f):
#                     gse = os.path.basename(gse_f).replace(".tsv", "")
#                     gse_passed_gq.append(gse)
#
#         with open(str(output), "w") as out:
#             for gse in gse_passed:
#                 print(gse, file=out)