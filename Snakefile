# TODO: create script to initiate files
import pandas as pd

configfile: "config.yaml"

#GSMS = [line.rstrip('\n') for line in open(config["gsm_list"])]
#GSE = ["GSE86552"]

GSMS = ["GSM1498954"]

SRR_TO_GSM = pd.read_csv("/home/boris/gq-rnaseq-pipeline/files/srr_to_gsm.tsv", sep="\t")

print(SRR_TO_GSM[SRR_TO_GSM.gsm == "GSM2305547"]["srr"].tolist())

rule all:
    input:
        expand(
            ["out/gsm/{gsm}.tsv"],
            gsm=GSMS
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
        sra = "out/sra/{srr}.sra",
        ref = config["refseq"]
    output:
        h5 = "out/kallisto/{srr}/abundance.h5",
        tsv = "out/kallisto/{srr}/abundance.tsv",
        json = "out/kallisto/{srr}/run_info.json"
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
    output: "out/gsm/{gsm}.tsv"
    log: "out/gsms/{gsm}.log"
    message: "Agregating GSM"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell: "Rscript ./scripts/srr_to_gsm.R {wildcards.gsm} \
                ./files/srr_to_gsm.tsv ./files/probes_to_genes.tsv"

# checkpoint check_gsm:
#     input:
#         lambda wildcards: expand(
#             "{gsm}",
#             gsm=GSMS[wildcards.gse]
#             )
#     output:
#         qc_filter_table = "{gse}_filtration.tsv"
#     message: "GSM QC checkpoint"
#
# def gsm_to_gse_input_fun(**wildcards):
#     path = checkpoints.check_gsm.get(**wildcards).output.qc_filter_table
#     #df = pd.read_csv(path, sep="\t")
#     #bams = df[[answer2bool(v) for v in df['qc_passed_sample']]]["File"]
#     gsms_passed_qc = []
#     return gsms_passed_qc
#
# rule gsm_to_gse:
#     input: gsm_to_gse_input_fun
#     output:
#             gse="{gse}.tsv"

