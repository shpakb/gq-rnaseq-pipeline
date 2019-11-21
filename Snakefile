# TODO: create script to initiate files
import pandas as pd
import os

configfile: "config_rn.yaml"

#GSES = [line.rstrip('\n') for line in open(config["gse_list"])]
GSES = ["GSE61179"]

SRR_TO_GSM = pd.read_csv(config["srr_to_gsm"], sep="\t")
GSM_TO_GSE = pd.read_csv(config["gsm_to_gse"], sep="\t")

rule all:
    input:
        expand(
            ["out/gses/{gse}.tsv"],
            gse=GSES
        )

rule sra_download:
    resources:
        load=25,
        mem_ram=1
    output: temp("out/sra/{srr}.sra")
    log:    "out/sra/{srr}.log"
    message: "Downloading {wildcards.srr}"
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/download_sra.sh {wildcards.srr} {output} >{log} 2>&1"

rule sra_kallisto_quant:
    input:
        sra="out/sra/{srr}.sra",
        ref=config["refseq"]
    output:
        h5="out/kallisto/{srr}/abundance.h5",
        tsv="out/kallisto/{srr}/abundance.tsv",
        json="out/kallisto/{srr}/run_info.json"
    log:
        "out/kallisto/{srr}/{srr}.log"

    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/quantify.sh {wildcards.srr} {input.sra} {input.ref} out/kallisto/{wildcards.srr} >{log} 2>&1"

rule srr_to_gsm:
    input:
        lambda wildcards: expand(
            "out/kallisto/{srr}/abundance.tsv",
            srr=SRR_TO_GSM[SRR_TO_GSM.gsm == wildcards.gsm]["srr"].tolist()
        )
    output: "out/gsms/{gsm}.tsv"
    log: "out/gsms/{gsm}.log"
    message: "Aggregating GSM"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/srr_to_gsm.R {wildcards.gsm}"
        " {config[probes_to_genes]} {config[srr_to_gsm]}"
        " out/kallisto out/gsms"

# Checkpoint
# If GSE is not empty appends GSE id to gse_passed_qc
checkpoint aggregate_filters:
    input:
        expand(
            "out/gses/{gse}.tsv",
            gse=GSES
        )
    output:
        "out/gse_passed_qc.tsv"
    run:
        gse_passed = []
        for gse_f in input:
            if os.path.getsize(gse_f):
                gse = os.path.basename(gse_f).replace(".tsv", "")
                gse_passed.append(gse)

        with open(str(output), "w") as out:
            for gse in gse_passed:
                print(gse, file=out)

# Aggregates GSMs to GSE
# If GSE fails QC outputs empty GSE
rule gsm_to_gse:
    input:
        lambda wildcards: expand(
            "out/{gsm}.tsv",
            gsm=GSM_TO_GSE[GSM_TO_GSE.gse == wildcards.gse]["gsm"].tolist()
        )
    output: "out/gses/{gse}.tsv"
    log: "out/gses/{gse}.log"
    message: "Aggregating GSE"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/gsm_to_gse.R {wildcards.gse}"
        " {config[probes_to_genes]} {config[gsm_to_gse]}"
        " out/gsms out/gses"
