import os

ALL_GSES=["GSE1", "GSE2"]
GSMS = {
    "GSE1": ["GSM10", "GSM11"],
    "GSE2": ["GSM20", "GSM21"],
    "GSE3": ["GSM3"],
    "GSE4": ["GSM3"],
}


rule all:
    input: "gse_passed_qc.tsv"

checkpoint aggregate_filters:
    input:
        expand(
            "gses/{gse}.tsv",
            gse=ALL_GSES
        )
    output: "gse_passed_qc.tsv"
    run:
        gse_passed = []
        for gse_f in input:
            if os.path.getsize(gse_f):
                gse = os.path.basename(gse_f).replace(".tsv", "")
                gse_passed.append(gse)

        with open(str(output), "w") as out:
            for gse in gse_passed:
                print(gse, file=out)


rule aggregate_gse:
    input:
        lambda wildcards: expand(
            "gsms/{gsm}.tsv",
            gsm=GSMS[wildcards.gse]
        )
    output:
        # filter = "{gse}_filter.tsv",
        gse = "gses/{gse}.tsv"
    run:
        if wildcards.gse == "GSE1":
            shell("touch {output}")
        else:
            shell("echo data > {output}")

rule quant:
    #input: expand sras...
    output: touch("gsms/{gsm}.tsv")