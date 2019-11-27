configfile: "configs/config_rn.yaml"

rule series_matrices_seq_download:
    '''
    Takes search output .txt, parses out all links and downloads them. 
    (wget -nc) Existing up to date files out dir doesn't reload. To update  database just remove complete flag file. 
    '''
    input:
        config["geo_search_result_chip"]
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        sm_dir=directory("out/series_matrices/seq/"),
        complete_flag="out/flags/series_matrices_seq_download"
    shell:
        "scripts/bash/download_sm.sh {input} {output.sm_dir} && "
        "touch {output.complete_flag}"

rule sm_seq_metadata:
    input:
        rules.series_matrices_seq_download.output.sm_dir
    output:
        gsm_table="out/data/metadata/seq/gsm.tsv",
        gse_table="out/data/metadata/seq/gse.tsv"
    shell:
        "python scripts/python/parse_sm_metadata.py {input} {output.gse_table} {output.gsm_table}"

rule sra_accession_table_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        temp("out/data/sra_accession_raw.tab")
    shell:
        "wget -O {output} {config[sra_accession_table]}"

rule get_srr_table:
    resources:
        writing_res=1
    input:
        rules.sra_accession_table_download.output
    output:
        srr_table="out/data/srr_gsm_spots.tsv"
    shell:
        "scripts/bash/clean_sra_accession_table.sh {input} {output.srr_table}"

rule prefilter_seq_gse:
    '''
    Takes seq SM metadata and prefiltering values from config (
    Outputs:
        1) gse_filtering.tsv
        2) gsm_filtering.tsv
        3) gse_prefiltered.list with list of GSE ids that has number of GSM
    between mim and max number of GSM after filtering out other types of experiments
    '''
    input:
        srr_table=rules.get_srr_table.output.srr_table,
        gse_table=rules.sm_seq_metadata.output.gse_table,
        gsm_table=rules.sm_seq_metadata.output.gsm_table
    output:
        gse_filtering_df="out/data/gse_filtering.tsv",
        gsm_filtering_df="out/data/gsm_filtering.tsv",
        gse_filtered_list="out/data/gse_filtered.list"
    shell:
        "Rscript scripts/R/filter_seq.R {input.srr_table} {input.gse_table} {input.gsm_table}"
rule sra_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output: temp("out/sra/{srr}.sra")
    log:    "out/sra/{srr}.log"
    message: "Downloading {wildcards.srr}"
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/bash/download_sra.sh {wildcards.srr} {output} > {log} 2>&1"

rule sra_fastqdump:
    resources:
        writing_res=1
    input:
        "out/sra/{srr}.sra"
    output:
        fastq_dir=temp(directory("out/fastq/{srr}")),
        complete_flag=temp("out/fastq/{srr}_complete")
    log:    "out/fastq/{srr}.log"
    message: "fastq-dump {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shell:
        "fastq-dump --outdir {output.fastq_dir} --split-3 {input} >{log} 2>&1 &&"
        "touch {output.complete_flag}"

rule fastq_kallisto:
    resources:
        mem_ram=4
    input:
        rules.sra_fastqdump.output.complete_flag,
        fastq_dir="out/fastq/{srr}"
    output:
        h5=protected("out/kallisto/{srr}/abundance.h5"),
        tsv=protected("out/kallisto/{srr}/abundance.tsv"),
        json=protected("out/kallisto/{srr}/run_info.json")
    log: "out/kallisto/{srr}/{srr}.log"
    message: "Kallisto: {wildcards.srr}"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/R/quantify.sh {wildcards.srr} {input.fastq_dir} "
        " {config[refseq]} out/kallisto/{wildcards.srr} >{log} 2>&1"

# rule srr_to_gsm:
#     resources:
#         mem_ram=1
#     input:
#         lambda wildcards: expand(
#             "out/kallisto/{srr}/abundance.tsv",
#             srr=SRR_GSM_DF[SRR_GSM_DF.gsm == wildcards.gsm]["srr"].tolist()
#         )
#     output: "out/gsms/{gsm}.tsv"
#     log: "out/gsms/{gsm}.log"
#     message: "Aggregating GSM {wildcards.gsm}"
#     shadow: "shallow"
#     conda: "envs/r_scripts.yaml"
#     shell:
#         "Rscript scripts/srr_to_gsm.R {wildcards.gsm}"
#         " {config[probes_to_genes]} {config[srr_to_gsm]}"
#         " out/kallisto out/gsms"

# checkpoint get_gsm_stats:
#     resources:
#         mem_ram=2
#     input:
#         gsm_files=lambda wildcards: expand(
#             "out/gsms/{gsm}.tsv",
#             gsm=GSM_GSE_DF[GSM_GSE_DF['gse'].isin(GSES)]["gsm"].tolist()
#         ),
#     output:
#         gsm_stats="out/stats/gsm_stats.tsv",
#         gse_filtered_df="out/stats/gse_gsm_filtered.tsv"
#     shell:
#         "Rscript scripts/get_gsm_stats.R {output.gsm_stats} {input.gsm_files} &&"
#         "   Rscript scripts/filtered_gse_list.R {output.gsm_stats}"
#         "       {output.gse_filtered_df} {config[gsm_to_gse]} {config[min_gse_gq]} {config[min_exp_genes]}"

rule gsm_to_gse:
    input:
        gse_filtered_df="out/stats/gse_gsm_filtered.tsv"
    output:
        gse="out/gses/{gse}.tsv"
    log: "out/gses/{gse}.log"
    message: "Aggregating GSE"
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/gsm_to_gse.R {output} out/gsms"
        " {config[ensamble_genesymbol_entrez]} "
        " {input.gse_filtered_df}"

# def get_filtered_gses(wildcards):
#     gse_filtered_path = checkpoints.get_gsm_stats.get(**wildcards).output.gse_filtered_df
#     gse_gsm_df = pd.read_csv(gse_filtered_path, sep="\t")
#     gses = list(set(gse_gsm_df["gse"].tolist()))
#     gse_files = ["out/gses/" + gse + ".tsv" for gse in gses]
#     return gse_files

# rule push_filtered_gses:
#     input:
#         get_filtered_gses
#     output:
#         flag="out/flag"
#     shell:
#         "touch {output.flag}"