#
#
# checkpoint select_gse:
#     input: "pass"
#     output: "pass"
#     run: pass
#
# rule wgcna_gse:
#     '''
#     Performs WGCNA on both chip and seq GSE. Outputs folder with two file for each
#     '''
#     input: "pass"
#     output: "pass"
#     run: pass
#
# checkpoint wgcna_stats:
#     '''
#         Parses WGCNA logs and gathers WGCNA clusterization stats to df.
#         STATUS- ok/failed+{reason}...
#         POWER, FIT, N_CLUSTERS
#     '''
#     input: "pass"
#     output: "out/data/wgcna/wgcna_stats.tsv"
#     run: pass
#
# def get_complete_wgcna:
#     '''
#     Returns list of out/wgcna/{gse}/modules.tsv
#     '''
#     return pass
#
# rule aggregate_gqdb:
#     '''
#     Aggregates
#     '''
#     input: "pass"
#     output: "out/data/wgcna/wgcna_stats.tsv"
#     run: pass
#
# rule get_eigengene_annotations:
#     '''
#     Outputs eigenegene heatmap annotations for dynamic heatmaps.
#     BUT HAS TO IMPLEMENT DYNAMIC HEATMAPS FIRST
#     '''
#     input: "pass"
#     output: "out/data/wgcna/wgcna_stats.tsv"
#     run: pass