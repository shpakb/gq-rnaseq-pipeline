# TODO: create script to initiate pipeline_wdir
configfile: "config.yaml"

GSES = [line.rstrip('\n') for line in open(config["gse_list"])]

rule all:
    input:
        expand(

                ["gses/{gse}.tsv"]
            ,
            gses=GSES
        )

rule sra_download:
    resources:
        load=25,
        mem_ram=1
    output: temp("sra/{srr}.sra")
    log:    "sra/{srr}.log"
    message: "Downloading {wildcards.srr}"
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/download_sra.sh {wildcards.srr} {output} >{log} 2>&1"


rule sra_kallisto_quant:
    input:
        sra = "sra/{srr}.sra",
        ref = config["refseq"]
    output:
        h5 = "kallisto/{srr}/abundance.h5",
        tsv = "kallisto/{srr}/abundance.tsv",
        json = "kallisto/{srr}/run_info.json"
    log:
        "kallisto/{srr}/{srr}.log"

    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/quantify.sh {wildcards.srr} {input.sra} {input.ref} kallisto/{wildcards.srr} >{log} 2>&1"

rule srr_to_gsm:
    input:
        expand(
            srr="kallisto/{srr}/abundance.tsv",
            srr=[]#TODO: write filter expression that takes GSM and outputs list of SRR
            )
    log: "gsms/{gsm}.log"
    message: "Agregating GSM"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell: "Rscript ./scripts/srr_to_gsm.R {wildcards.gsm} \
                ./files/srr_to_gsm.tsv ./files/probes_to_genes.tsv"

checkpoint check_gsm:
    input:
        lambda wildcards: expand(
            "{gsm}",
            gsm=GSMS[wildcards.gse]
            )
    output:
        qc_filter_table = "{gse}_filtration.tsv"
    # filtered gsm

def gsm_to_gse_input_fun(**wildcards):
    path = checkpoints.check_gsm.get(**wildcards).output.qc_filter_table
    #df = pd.read_csv(path, sep="\t")
    #bams = df[[answer2bool(v) for v in df['qc_passed_sample']]]["File"]
    gsms_passed_qc = []
    return gsms_passed_qc

rule gsm_to_gse:
    input: gsm_to_gse_input_fun
    output:
            gse="{gse}.tsv"


# Find rat GSM made of few SRR