
# checkpoint filter_chip_gse:
#     '''
#     Takes chip sm metadata
#     Outputs {lab}_filtered.list with list of GSE ids that:
#         1) Has exp df
#         2) GPL in list of accepted df platforms be annotated with 3col
#     '''
#     output: "out"
#
#
# rule series_matrices_chip_download:
#     '''
#     Parses search output .txt for all links and downloads them.
#     wget -nc; Downloaded files in folder doesn't reload. To update
#     database just remove complete flag file.
#     '''
#     input:
#         lambda wildcards: config[geo_search_result_chip]
#     resources:
#         download_res=1,
#         writing_res=1,
#         mem_ram=2
#     output:
#         sm_dir=directory("out/series_matrices/chip"),
#         complete_flag="out/flags/series_matrices_chip_download_complete"
#     shell:
#         "scripts/bash/download_sm.sh {output.sm_dir} && "
#         "touch {output.complete_flag}"