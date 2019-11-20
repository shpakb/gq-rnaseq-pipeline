configfile: "config.yaml"

SRAS = [line.rstrip('\n') for line in open(config["srr_list"])]

rule all:
    input:
        expand(
            [
                "kallisto/{srr}/abundance.h5",
                "kallisto/{srr}/abundance.tsv",
                "kallisto/{srr}/run_info.json"
            ],
            srr=SRAS
        )

rule sra_download:
    resources: load=25
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