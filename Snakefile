
include: "raw_data_seq.smk"
include: "raw_data_chip.smk"
include: "gds_check.smk"
include: "gq_processing.smk"
include: "pca_processing.smk"


rule all:
    input:
        rules.series_matrices_seq_download.output.sm_dir,
        rules.series_matrices_seq_download.output.complete_flag,
        rules.sm_seq_metadata.output.gsm_table,
        rules.sm_seq_metadata.output.gse_table,
        rules.get_srr_gsm_spots_table.output
