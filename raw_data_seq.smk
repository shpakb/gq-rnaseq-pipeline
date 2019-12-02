# Issues: can't pass files by input functions from checkpoints in runtime. But still works during graph calculation.
# 1) Cant pass temp(directory("out/fasqdump/")) as reference to from rules variables. Had to make workaround with
# completion flags

# Add list with priority GSE that bypasses both checkpoints
# (except for the filtering on right organism and right type of samples)

# TODO: move all QC logic to input functions(if possible)


# rule sra_download:
#     resources:
#         download_res=1,
#         writing_res=1,
#         mem_ram=2
#     priority: 3
#     output: temp("out/{organism}/sra/{srr}.sra")
#     log:    "logs/{organism}/sra_download/{srr}.log"
#     message: "Downloading {wildcards.srr}"
#     shadow: "shallow"
#     conda: "envs/quantify.yaml"
#     shell:
#         "scripts/bash/download_sra.sh {wildcards.srr} {output}"
#         " > {log} 2>&1"
#
# rule sra_fastqdump:
#     resources:
#         writing_res=1
#     input:
#         "out/sra/{srr}.sra"
#     output:
#         fastq_dir=temp(directory("out/{organism}/fastq/{srr}")),
#         complete_flag=temp("out/{organism}/fastq/{srr}_complete")
#     log:    "logs/{organism}/sra_fastqdump/{srr}.log"
#     message: "fastq-dump {wildcards.srr}"
#     conda: "envs/quantify.yaml"
#     shell:
#         "fastq-dump --outdir {output.fastq_dir} --split-3 {input}"
#         " > {log} 2>&1 &&"
#         " touch {output.complete_flag}"
#
# rule fastq_kallisto:
#     resources:
#         mem_ram=lambda wildcards: {"hs": 10, "mm": 8, "rn": 5}[{wildcards.organism}]
#     priority: 2
#     input:
#         rules.sra_fastqdump.output.complete_flag,
#         fastq_dir="out/{organism}/seq/fastq/{srr}",
#         refseq="input/{organism}/seq/refseq_kallisto"
#     output:
#         h5=protected("out/{organism}/seq/kallisto/{srr}/abundance.h5"),
#         tsv=protected("out/{organism}/seq/kallisto/{srr}/abundance.tsv"),
#         json=protected("out/{organism}/seq/{srr}/run_info.json")
#     log: "logs/{organism}/seq/fastq_kallisto/{srr}.log"
#     message: "Kallisto: {wildcards.srr}"
#     conda: "envs/quantify.yaml"
#     shadow: "shallow"
#     shell:
#         "scripts/bash/quantify.sh {wildcards.srr} {input.fastq_dir} {config[refseq]} out/kallisto/{wildcards.srr}"
#         " > {log} 2>&1"
#
# # For some reason input function returns proper output during graph calculation but
# # returns first output variable from checkpoint during rule execution
# # This does guaranties that all necessary files for the task are there but doesn't allow to pass them as an input to
# # script.
#
# def get_srr_files(wildcards):
#     srr_df_file = checkpoints.prequant_filter.get(**wildcards).output.srr_gsm_df
#     srr_df = pd.read_csv(srr_df_file, sep="\t")
#     srr_list = srr_df[srr_df['GSM']==wildcards.gsm]["SRR"].tolist()
#     srr_files = expand("out/{wildcards.organism}/seq/kallisto/{srr}/abundance.tsv", srr=srr_list)
#     return srr_files
#
# rule srr_to_gsm:
#     resources:
#         mem_ram=1
#     input:
#         srr_files=get_srr_files,
#         srr_df="out/{organism}/seq/data/filtering/prequant/srr_gsm.tsv"
#     output:
#         gsm_file="out/{organism}/seq/gsms/{gsm}.tsv"
#     log: "logs/{organism}/srr_to_gsm/{gsm}.log"
#     message: "Aggregating {wildcards.gsm}"
#     shadow: "shallow"
#     conda: "envs/r_scripts.yaml"
#     shell:
#         "Rscript scripts/R/srr_to_gsm.R {output.gsm_file} {config[probes_to_genes]} {input.srr_df}"
#         " > {log} 2>&1"
#
# def get_prequant_filtered_gsm_files(wildcards):
#     filtered_gsm_list = checkpoints.prequant_filter.get(**wildcards).output.passing_gsm_list
#     gsm_list = [line.rstrip('\n') for line in open(filtered_gsm_list)]
#     gsm_files=expand("out/{wildcards.organism}/seq/gsms/{gsm}.tsv",gsm=gsm_list)
#     return gsm_files
#
# def get_gsm_gse_df(wildcards):
#     return checkpoints.prequant_filter.get(**wildcards).output.gsm_gse_df
#
# checkpoint postquant_filter:
#     input:
#         gsm_files=get_prequant_filtered_gsm_files,
#         gsm_gse_df=get_gsm_gse_df
#     output:
#         gsm_stats_df="out/{oragnism}/seq/data/filtering/postquant/gsm_stats.tsv",
#         gsm_gse_df="out/{oragnism}/seq/data/filtering/postquant/gsm_gse.tsv",
#         passing_gse_list="out/{oragnism}/seq/data/filtering/postquant/passing_gse.list"
#     message: "Post quantification filtering."
#     log: "logs/{oragnism}/postquant_filter.log"
#     conda: "envs/r_scripts.yaml"
#     shell:
#         "Rscript scripts/R/postquant_filter.R {config[quant_min_gsm]} {config[min_exp_genes]} {input.gsm_gse_df}"
#         " {output.gsm_stats_df} {output.gsm_gse_df} {output.passing_gse_list} {input.gsm_files}"
#         " > {log} 2>&1"
#
# def get_postquant_gsms_for_gse(wildcards):
#     gsm_gse_df_file = checkpoints.postquant_filter.get(**wildcards).output.gsm_gse_df
#     gsm_gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
#     gsms = gsm_gse_df[gsm_gse_df['GSE']==wildcards.gse]["GSM"].tolist()
#     gsm_files = expand("out/{wildcards.organism}/seq/gsms/{gsm}.tsv", gsm=gsms)
#     return gsm_files
#
# rule gsm_to_gse:
#     input:
#         gsm_files=get_postquant_gsms_for_gse,
#         gse_df="out/{organism}/seq/data/filtering/postquant/gsm_gse.tsv",
#         gene_mapping="input/{organism}/ensamble_genesymbol_entrez.tsv}"
#     output:
#         gse="out/{organism}/seq/gses/{gse}.tsv"
#     log: "logs/{organism}/gsm_to_gse/{gse}.log"
#     message: "Aggregating {wildcards.gse}"
#     shadow: "shallow"
#     conda: "envs/r_scripts.yaml"
#     shell:
#         "Rscript scripts/R/gsm_to_gse.R {output.gse} out/{wildcards.organism}/seq/gsms {input.gene_mapping} "
#         " {input.gse_df}"
#         " > {log} 2>&1"
#
# def get_postquant_passing_gse(wildcards):
#     filtered_gse_list = checkpoints.postquant_filter.get(**wildcards).output.passing_gse_list
#     gses = [line.rstrip('\n') for line in open(filtered_gse_list)]
#     gse_files = expand("out/{wildcards.organism}/seq/gses/{gsm}.tsv", gsm=gses)
#     return gse_files
#
# rule push_seq_gses:
#     input:
#         get_postquant_passing_gse
#     output:
#         flag="out/flags/{organism}_seq_complete"
#     shell:
#         "touch {output.flag}"
