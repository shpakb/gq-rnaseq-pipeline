# TODO: run Sashas dataset quantification

include: "raw_data_seq.smk"
include: "raw_data_chip.smk"
include: "gds_check.smk"
include: "gq_processing.smk"
include: "pca_processing.smk"


rule all:
    input:
        rules.push_filtered_gses.output.flag