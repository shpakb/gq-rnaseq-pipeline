configfile: "config.yaml"

# TODO: PCA query on rna-seq
# TODO: aggregate metadata on kallisto runs in df and choose version of callisto used for most of the samples.
# TODO: gz all tsv files in pipeline


import pandas as pd
from os import listdir
from os.path import isfile, join
import re

glob_srr_df = {}

rule all:
    input:
        expand("out/mm/seq/gsms/{gsm}.tsv", gsm=['GSM2373319',
                                                 'GSM2373313',
                                                 'GSM2373184',
                                                 'GSM2373480',
                                                 'GSM2373336',
                                                 'GSM2373358',
                                                 'GSM2373454',
                                                 'GSM2373362',
                                                 'GSM2373427',
                                                 'GSM2373195',
                                                 'GSM2373415',
                                                 'GSM2373253',
                                                 'GSM2373227',
                                                 'GSM2373346',
                                                 'GSM2373242',
                                                 'GSM2373177',
                                                 'GSM2373250',
                                                 'GSM2373429',
                                                 'GSM2373243',
                                                 'GSM2373378',
                                                 'GSM2373393',
                                                 'GSM2373315',
                                                 'GSM2373203',
                                                 'GSM2373390',
                                                 'GSM2373159',
                                                 'GSM2373145',
                                                 'GSM2373351',
                                                 'GSM2373370',
                                                 'GSM2373214',
                                                 'GSM2373174',
                                                 'GSM2373171',
                                                 'GSM2373361',
                                                 'GSM2373410',
                                                 'GSM2373123',
                                                 'GSM2373477',
                                                 'GSM2373489',
                                                 'GSM2373316',
                                                 'GSM2373169',
                                                 'GSM2373354',
                                                 'GSM2373292',
                                                 'GSM2373265',
                                                 'GSM2373231',
                                                 'GSM2373121',
                                                 'GSM2373364',
                                                 'GSM2373470',
                                                 'GSM2373274',
                                                 'GSM2373438',
                                                 'GSM2373269',
                                                 'GSM2373381',
                                                 'GSM2373468',
                                                 'GSM2373280',
                                                 'GSM2373262',
                                                 'GSM2373314',
                                                 'GSM2373251',
                                                 'GSM2373192',
                                                 'GSM2373191',
                                                 'GSM2373458',
                                                 'GSM2373447',
                                                 'GSM2373388',
                                                 'GSM2373295',
                                                 'GSM2373118',
                                                 'GSM2373401',
                                                 'GSM2373120',
                                                 'GSM2373366',
                                                 'GSM2373395',
                                                 'GSM2373368',
                                                 'GSM2373306',
                                                 'GSM2373113',
                                                 'GSM2373266',
                                                 'GSM2373329',
                                                 'GSM2373143',
                                                 'GSM2373455',
                                                 'GSM2373148',
                                                 'GSM2373387',
                                                 'GSM2373344',
                                                 'GSM2373365',
                                                 'GSM2373352',
                                                 'GSM2373487',
                                                 'GSM2373432',
                                                 'GSM2373416',
                                                 'GSM2373130',
                                                 'GSM2373288',
                                                 'GSM2373199',
                                                 'GSM2373116',
                                                 'GSM2373144',
                                                 'GSM2373210',
                                                 'GSM2373466',
                                                 'GSM2373252',
                                                 'GSM2373367',
                                                 'GSM2373249',
                                                 'GSM2373363',
                                                 'GSM2373296',
                                                 'GSM2373371',
                                                 'GSM2373208',
                                                 'GSM2373404',
                                                 'GSM2373287',
                                                 'GSM2373347',
                                                 'GSM2373150',
                                                 'GSM2373189',
                                                 'GSM2373348',
                                                 'GSM2373115',
                                                 'GSM2373142',
                                                 'GSM2373465',
                                                 'GSM2373211',
                                                 'GSM2373452',
                                                 'GSM2373239',
                                                 'GSM2373225',
                                                 'GSM2373493',
                                                 'GSM2373476',
                                                 'GSM2373275',
                                                 'GSM2373273',
                                                 'GSM2373428',
                                                 'GSM2373256',
                                                 'GSM2373293',
                                                 'GSM2373435',
                                                 'GSM2373212',
                                                 'GSM2373406',
                                                 'GSM2373125',
                                                 'GSM2373442',
                                                 'GSM2373175',
                                                 'GSM2373298',
                                                 'GSM2373282',
                                                 'GSM2373384',
                                                 'GSM2373270',
                                                 'GSM2373186',
                                                 'GSM2373309',
                                                 'GSM2373158',
                                                 'GSM2373369',
                                                 'GSM2373457',
                                                 'GSM2373236',
                                                 'GSM2373131',
                                                 'GSM2373247',
                                                 'GSM2373321',
                                                 'GSM2373338',
                                                 'GSM2373182',
                                                 'GSM2373173',
                                                 'GSM2373276',
                                                 'GSM2373446',
                                                 'GSM2373217',
                                                 'GSM2373155',
                                                 'GSM2373268',
                                                 'GSM2373301',
                                                 'GSM2373467',
                                                 'GSM2373204',
                                                 'GSM2373448',
                                                 'GSM2373258',
                                                 'GSM2373383',
                                                 'GSM2373323',
                                                 'GSM2373190',
                                                 'GSM2373303',
                                                 'GSM2373198',
                                                 'GSM2373202',
                                                 'GSM2373237',
                                                 'GSM2373311',
                                                 'GSM2373420',
                                                 'GSM2373343',
                                                 'GSM2373193',
                                                 'GSM2373445',
                                                 'GSM2373284',
                                                 'GSM2373219',
                                                 'GSM2373312',
                                                 'GSM2373324',
                                                 'GSM2373200',
                                                 'GSM2373430',
                                                 'GSM2373471',
                                                 'GSM2373494',
                                                 'GSM2373286',
                                                 'GSM2373134',
                                                 'GSM2373290',
                                                 'GSM2373320',
                                                 'GSM2373439',
                                                 'GSM2373221',
                                                 'GSM2373122',
                                                 'GSM2373126',
                                                 'GSM2373168',
                                                 'GSM2373359',
                                                 'GSM2373443',
                                                 'GSM2373441',
                                                 'GSM2373291',
                                                 'GSM2373241',
                                                 'GSM2373424',
                                                 'GSM2373353',
                                                 'GSM2373374',
                                                 'GSM2373117',
                                                 'GSM2373473',
                                                 'GSM2373350',
                                                 'GSM2373409',
                                                 'GSM2373399',
                                                 'GSM2373135',
                                                 'GSM2373220',
                                                 'GSM2373356',
                                                 'GSM2373235',
                                                 'GSM2373147',
                                                 'GSM2373418',
                                                 'GSM2373449',
                                                 'GSM2373222',
                                                 'GSM2373440',
                                                 'GSM2363332',
                                                 'GSM2373185',
                                                 'GSM2373163',
                                                 'GSM2373178',
                                                 'GSM2373332',
                                                 'GSM2373349',
                                                 'GSM2373172',
                                                 'GSM2373322',
                                                 'GSM2373170',
                                                 'GSM2373488',
                                                 'GSM2373183',
                                                 'GSM2373133',
                                                 'GSM2373464',
                                                 'GSM2373394',
                                                 'GSM2373188',
                                                 'GSM2373310',
                                                 'GSM2373339',
                                                 'GSM2373167',
                                                 'GSM2373114',
                                                 'GSM2373149',
                                                 'GSM2373146',
                                                 'GSM2373376',
                                                 'GSM2373472',
                                                 'GSM2373299',
                                                 'GSM2373257',
                                                 'GSM2373304',
                                                 'GSM2373283',
                                                 'GSM2373492',
                                                 'GSM2373485',
                                                 'GSM2373140',
                                                 'GSM2373111',
                                                 'GSM2373233',
                                                 'GSM2363333',
                                                 'GSM2373165',
                                                 'GSM2373128',
                                                 'GSM2373330',
                                                 'GSM2373461',
                                                 'GSM2373261',
                                                 'GSM2373490',
                                                 'GSM2373207',
                                                 'GSM2373218',
                                                 'GSM2373132',
                                                 'GSM2373307',
                                                 'GSM2373160',
                                                 'GSM2373232',
                                                 'GSM2373414',
                                                 'GSM2373259',
                                                 'GSM2373345',
                                                 'GSM2373176',
                                                 'GSM2373271',
                                                 'GSM2373196',
                                                 'GSM2373423',
                                                 'GSM2373360',
                                                 'GSM2373408',
                                                 'GSM2373154',
                                                 'GSM2373281',
                                                 'GSM2373333',
                                                 'GSM2373436',
                                                 'GSM2373469',
                                                 'GSM2373412',
                                                 'GSM2373335',
                                                 'GSM2373450',
                                                 'GSM2373161',
                                                 'GSM2373417',
                                                 'GSM2373318',
                                                 'GSM2373228',
                                                 'GSM2373342',
                                                 'GSM2373277',
                                                 'GSM2373124',
                                                 'GSM2373411',
                                                 'GSM2373327',
                                                 'GSM2373434',
                                                 'GSM2373285',
                                                 'GSM2373157',
                                                 'GSM2373437',
                                                 'GSM2373397',
                                                 'GSM2373459',
                                                 'GSM2373255',
                                                 'GSM2363334',
                                                 'GSM2373341',
                                                 'GSM2373396',
                                                 'GSM2373391',
                                                 'GSM2373141',
                                                 'GSM2373334',
                                                 'GSM2373479',
                                                 'GSM2373166',
                                                 'GSM2373272',
                                                 'GSM2373156',
                                                 'GSM2373264',
                                                 'GSM2373377',
                                                 'GSM2373413',
                                                 'GSM2373385',
                                                 'GSM2373453',
                                                 'GSM2373357',
                                                 'GSM2373398',
                                                 'GSM2373294',
                                                 'GSM2373254',
                                                 'GSM2373151',
                                                 'GSM2373201',
                                                 'GSM2373317',
                                                 'GSM2373300',
                                                 'GSM2373422',
                                                 'GSM2373267',
                                                 'GSM2373245',
                                                 'GSM2373389',
                                                 'GSM2373136',
                                                 'GSM2373164',
                                                 'GSM2373337',
                                                 'GSM2373179',
                                                 'GSM2373240',
                                                 'GSM2373153',
                                                 'GSM2373405',
                                                 'GSM2373139',
                                                 'GSM2373386',
                                                 'GSM2373379',
                                                 'GSM2373112',
                                                 'GSM2373229',
                                                 'GSM2373331',
                                                 'GSM2373481',
                                                 'GSM2373205',
                                                 'GSM2373244',
                                                 'GSM2373380',
                                                 'GSM2373289',
                                                 'GSM2373328',
                                                 'GSM2373403',
                                                 'GSM2373326',
                                                 'GSM2373340',
                                                 'GSM2373475',
                                                 'GSM2373180',
                                                 'GSM2373419',
                                                 'GSM2373460',
                                                 'GSM2373491',
                                                 'GSM2373483',
                                                 'GSM2373478',
                                                 'GSM2373425',
                                                 'GSM2373382',
                                                 'GSM2373238',
                                                 'GSM2373119',
                                                 'GSM2373215',
                                                 'GSM2373392',
                                                 'GSM2373431',
                                                 'GSM2373260',
                                                 'GSM2373197',
                                                 'GSM2373355',
                                                 'GSM2373187',
                                                 'GSM2373216',
                                                 'GSM2373482',
                                                 'GSM2373224',
                                                 'GSM2373213',
                                                 'GSM2373451',
                                                 'GSM2373407',
                                                 'GSM2373325',
                                                 'GSM2373486',
                                                 'GSM2373234',
                                                 'GSM2373226',
                                                 'GSM2373302',
                                                 'GSM2373181',
                                                 'GSM2373373',
                                                 'GSM2373421',
                                                 'GSM2373444',
                                                 'GSM2373372',
                                                 'GSM2373375',
                                                 'GSM2373206',
                                                 'GSM2373474',
                                                 'GSM2373433',
                                                 'GSM2373278',
                                                 'GSM2373162',
                                                 'GSM2373297',
                                                 'GSM2373400',
                                                 'GSM2373279',
                                                 'GSM2373152',
                                                 'GSM2373209',
                                                 'GSM2373230',
                                                 'GSM2373462',
                                                 'GSM2373127',
                                                 'GSM2373463',
                                                 'GSM2373484',
                                                 'GSM2373137',
                                                 'GSM2373138',
                                                 'GSM2373402',
                                                 'GSM2373223',
                                                 'GSM2373426',
                                                 'GSM2373194',
                                                 'GSM2373308',
                                                 'GSM2373246',
                                                 'GSM2373129',
                                                 'GSM2373263',
                                                 'GSM2373248',
                                                 'GSM2373305',
                                                 'GSM2373456'])
        #"out/mm/seq/postquant_filter/gsm_stats.tsv"
        # expand("flags/{organism}/{platform}/pca/{n_genes}_{scale}/flag",
        #         organism=['mm'],
        #         platform=['chip'],
        #         n_genes=config['pca_n_genes'],
        #         scale=config['pca_scale'])
        # expand("out/{organism}/{platform}/pca_fgsea/"
        #        "{max_genes}_{scale}_{max_comp}_{var_threshold}/"
        #        "prepared/{geneset_name}.tsv",
        #     organism=['mm'],
        #     platform=['chip'],
        #     max_genes=config['pca_n_genes'],
        #     scale=config['pca_scale'],
        #     max_comp="10",
        #     var_threshold="0.02",
        #     geneset_name=['HALLMARK_HYPOXIA']# 'HALLMARK_PANCREAS_BETA_CELLS', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING']
        #                   #'HALLMARK_SPERMATOGENESIS', 'HALLMARK_FATTY_ACID_METABOLISM', 'HALLMARK_BILE_ACID_METABOLISM',
        #                   # 'HALLMARK_P53_PATHWAY', 'HALLMARK_MYOGENESIS', 'HALLMARK_PROTEIN_SECRETION',
        #                   # 'HALLMARK_UV_RESPONSE_DN', 'HALLMARK_ANGIOGENESIS', 'HALLMARK_NOTCH_SIGNALING',
        #                   # 'HALLMARK_MYC_TARGETS_V2', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB', 'HALLMARK_KRAS_SIGNALING_DN',
        #                   # 'HALLMARK_HEDGEHOG_SIGNALING', 'HALLMARK_APICAL_SURFACE', 'HALLMARK_MYC_TARGETS_V1',
        #                   # 'HALLMARK_ALLOGRAFT_REJECTION', 'HALLMARK_CHOLESTEROL_HOMEOSTASIS',
        #                   # 'HALLMARK_ANDROGEN_RESPONSE', 'HALLMARK_E2F_TARGETS', 'HALLMARK_GLYCOLYSIS',
        #                   # 'HALLMARK_DNA_REPAIR', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
        #                   # 'HALLMARK_IL6_JAK_STAT3_SIGNALING', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
        #                   # 'HALLMARK_UNFOLDED_PROTEIN_RESPONSE', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
        #                   # 'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_UV_RESPONSE_UP',
        #                   # 'HALLMARK_WNT_BETA_CATENIN_SIGNALING', 'HALLMARK_INTERFERON_ALPHA_RESPONSE',
        #                   # 'HALLMARK_G2M_CHECKPOINT', 'HALLMARK_IL2_STAT5_SIGNALING', 'HALLMARK_APOPTOSIS',
        #                   # 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_ESTROGEN_RESPONSE_LATE',
        #                   # 'HALLMARK_COAGULATION', 'HALLMARK_XENOBIOTIC_METABOLISM', 'HALLMARK_COMPLEMENT',
        #                   # 'HALLMARK_ADIPOGENESIS', 'HALLMARK_TGF_BETA_SIGNALING', 'HALLMARK_MITOTIC_SPINDLE',
        #                   # 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_APICAL_JUNCTION', 'HALLMARK_KRAS_SIGNALING_UP',
        #                   # 'HALLMARK_PEROXISOME', 'HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_HEME_METABOLISM']
        # )
rule sm_download:
    '''
    Takes search output .txt, parses out all links and downloads them. (wget -nc) Existing up to date files out dir
    doesn't reload. Writes completion flag in the end.
    '''
    resources:
        time=60*24*2,
        download_res=1,
        writing_res=1,
        priority=5
    input:
        "input/{organism}/{platform}/gds_search_result.txt"
    output:
        "flags/{organism}/{platform}/sm_download.flag"
    params:
        sm_download_dir="out/{organism}/{platform}/series_matrices"
    message: "Downloading series matrices for {wildcards.organism} {wildcards.platform}..."
    log: "logs/{organism}/{platform}/sm_download.log"
    shell:
        "scripts/bash/download_sm.sh {input} {params.sm_download_dir} "
        " > {log} 2>&1 &&"
        " touch {output}"

checkpoint extract_sm_metadata:
    input:
        rules.sm_download.output
    output:
        gsm_df="out/{organism}/{platform}/sm_metadata/gsm.tsv",
        gse_df="out/{organism}/{platform}/sm_metadata/gse.tsv"
    params:
        sm_download_dir=rules.sm_download.params.sm_download_dir
    message: "Aggregating metadata from series matrices {wildcards.organism} {wildcards.platform} ..."
    log: "logs/{organism}/{platform}/sm_metadata.log"
    shell:
        "python scripts/python/parse_sm_metadata.py {params.sm_download_dir} {output.gse_df} {output.gsm_df}"
        " > {log} 2>&1"

#############################################RNASEQ#####################################################################
rule sra_accession_df_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    output:
        temp("out/data/sra_accession_raw.tab")
    message: "Downloading SRA accession table..."
    log: "logs/sra_accession_df_download.log"
    shell:
        "wget -O {output} {config[sra_accession_df]}"
        " > {log} 2>&1"

rule get_srr_df:
    resources:
        writing_res=1
    input:
        rules.sra_accession_df_download.output
    output:
        "out/data/srr_gsm_spots.tsv"
    message: "Cleaning SRR df..."
    log: "logs/get_srr_df.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "scripts/bash/clean_sra_accession_df.sh {input} {output}"
        " > {log} 2>&1 &&"
        " Rscript scripts/R/clean_sra_accession_df.R {output}"
        " >> {log} 2>&1"

checkpoint prequant_filter:
    '''
    Outputs QC tables and lists of passing GSE and GSM. Takes one extra argument: priority_gse_list with GSEs that get
    to downstream quantification regardless of QC. Keep in mind that GSM for priority GSE gets filtered by sample
    metadata(GSM can be microarray, wrong organism etc.) 
    Temporary added priority_only flag and script logic for test runs and to quantify things while pipeline under 
    development.
    '''
    input:
        srr_df=rules.get_srr_df.output,
        gse_df="out/{organism}/seq/sm_metadata/gse.tsv",
        gsm_df="out/{organism}/seq/sm_metadata/gsm.tsv",
        gpl_df="input/gpl.tsv",
        priority_gse_list="input/{organism}/seq/priority_gse.list"
    params:
        min_spots=lambda wildcards: config["min_spots"][wildcards.organism],
        max_spots=lambda wildcards: config["max_spots"][wildcards.organism],
        organism=lambda wildcards: config["organism"][wildcards.organism]
    output:
        gsm_gse_df="out/{organism}/seq/prequant_filter/gsm_gse.tsv",
        gsm_filtering_df="out/{organism}/seq/prequant_filter/gsm_filtering.tsv",
        passing_gsm_list="out/{organism}/seq/prequant_filter/passing_gsm.list",
        srr_gsm_df="out/{organism}/seq/prequant_filter/srr_gsm.tsv"
    message: "Pre-quantification filtering for {wildcards.organism}..."
    log: "logs/{organism}/seq/prequant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/prequant_filter.R {input.gse_df} {input.gsm_df} {input.gpl_df} {input.srr_df}"
        " {params.organism} {params.min_spots} {params.max_spots} {config[quant_min_gsm]}"
        " {config[quant_max_gsm]} {output.gsm_filtering_df} {output.passing_gsm_list} {output.srr_gsm_df}"
        " {output.gsm_gse_df} {input.priority_gse_list} {config[priority_only_f]}"
        " > {log} 2>&1"

rule sra_download:
    resources:
        download_res=1,
        writing_res=1,
        mem_ram=2
    priority: 3
    output: temp("out/{organism}/seq/sra/{srr}.sra")
    log:    "logs/{organism}/seq/sra_download/{srr}.log"
    message: "Downloading {wildcards.srr} ({wildcards.organism})..."
    shadow: "shallow"
    conda: "envs/quantify.yaml"
    shell:
        "scripts/bash/download_sra.sh {wildcards.srr} {output}"
        " > {log} 2>&1"

rule sra_fastqdump:
    resources:
        writing_res=1
    input:
        rules.sra_download.output
    output:
        fastq_dir=temp(directory("out/{organism}/seq/fastq/{srr}")),
    log:    "logs/{organism}/seq/sra_fastqdump/{srr}.log"
    message: "fastq-dump {wildcards.srr} ({wildcards.organism})"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "fastq-dump --outdir {output.fastq_dir} --split-3 {input}"
        " > {log} 2>&1"

rule fastq_kallisto:
    resources:
        mem_ram=lambda wildcards: config["quant_mem_ram"][wildcards.organism]
    priority: 2
    input:
        fastq_dir=rules.sra_fastqdump.output,
        refseq="input/{organism}/seq/refseq_kallisto"
    output:
        # before uncommenting h5 remove files of older versions as described in README. See NOTES on h5
        #h5=protected("out/{organism}/seq/kallisto/{srr}/abundance.h5"),
        tsv=protected("out/{organism}/seq/kallisto/{srr}/abundance.tsv"),
        json=protected("out/{organism}/seq/kallisto/{srr}/run_info.json")
    log: "logs/{organism}/seq/fastq_kallisto/{srr}.log"
    message: "Kallisto: {wildcards.srr} ({wildcards.organism})"
    conda: "envs/quantify.yaml"
    shadow: "shallow"
    shell:
        "scripts/bash/quantify.sh {wildcards.srr} {input.fastq_dir} {input.refseq}"
        " out/{wildcards.organism}/seq/kallisto/{wildcards.srr}"# output dir 
        " {config[n_bootstrap]} " 
        " > {log} 2>&1"

def get_srr_for_gsm(wildcards):
    print(f"Getting SRR for {wildcards.gsm}")
    global glob_srr_df

    # for speed. Otherwise it takes forever
    if wildcards.organism not in glob_srr_df:
        print("Saving srr df to global var")
        srr_df_file = f"out/{wildcards.organism}/seq/prequant_filter/srr_gsm.tsv"
        srr_df = pd.read_csv(srr_df_file, sep="\t")
        glob_srr_df[wildcards.organism] = srr_df

    srr_df = glob_srr_df[wildcards.organism]
    srr_list = srr_df[srr_df['GSM']==wildcards.gsm]["SRR"].tolist()
    srr_list = list(set(srr_list)) # removing possible duplicates

    #srr_files = [checkpoints.fastq_kallisto.get(organism=wildcards.organism, srr=srr).output.tsv for srr in srr_list]

    srr_files = \
        expand("out/{organism}/seq/kallisto/{srr}/abundance.tsv",
        srr=srr_list,
        organism=wildcards.organism)

    return srr_files

rule srr_to_gsm:
    resources:
        mem_ram=1
    input:
        srr_files=get_srr_for_gsm,
        transcript_gene="input/{organism}/seq/transcript_gene.tsv"
    output:
        "out/{organism}/seq/gsms/{gsm}.tsv"
    log: "logs/{organism}/seq/srr_to_gsm/{gsm}.log"
    message: "Aggregating {wildcards.gsm} ({wildcards.organism})..."
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/srr_to_gsm.R {output} {input.transcript_gene} {input.srr_files}"
        " > {log} 2>&1"

def prequant_filtered_gsm_files(wildcards):
    print("Getting prequant filtered GSM files")
    filtered_gsm_list = checkpoints.prequant_filter.get(**wildcards).output.passing_gsm_list
    gsm_list = [line.rstrip('\n') for line in open(filtered_gsm_list)]
    #gsm_files = [checkpoints.srr_to_gsm.get(organism=wildcards.organism, gsm=gsm).output.tsv for gsm in gsm_list]
    gsm_files=\
        expand("out/{organism}/seq/gsms/{gsm}.tsv",
            gsm=gsm_list,
            organism=wildcards.organism)
    return gsm_files

checkpoint postquant_filter:
    input:
        gsm_files=prequant_filtered_gsm_files,
        gsm_gse_df=lambda wildcards: checkpoints.prequant_filter.get(**wildcards).output.gsm_gse_df
    output:
        gsm_stats_df="out/{organism}/seq/postquant_filter/gsm_stats.tsv",
        gsm_gse_df="out/{organism}/seq/postquant_filter/gsm_gse.tsv",
        passing_gse_list="out/{organism}/seq/postquant_filter/passing_gse.list"
    message: "Post quantification filtering ({wildcards.organism}) ..."
    params:
        min_exp_genes=lambda wildcards: config["min_exp_genes"][wildcards.organism]
    log: "logs/{organism}/seq/postquant_filter.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/postquant_filter.R {config[quant_min_gsm]} {params.min_exp_genes} {input.gsm_gse_df}"
        " {output.gsm_stats_df} {output.gsm_gse_df} {output.passing_gse_list} {input.gsm_files}"
        " > {log} 2>&1"

def postquant_gsms_for_gse(wildcards):
    print(f"Getting postquant filtered GSM files for {wildcards.gse}")
    gsm_gse_df_file = f"out/{wildcards.organism}/seq/postquant_filter/gsm_gse.tsv"
    gsm_gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
    gsms = gsm_gse_df[gsm_gse_df['GSE']==wildcards.gse]["GSM"].tolist()
    gsm_files = \
        expand("out/{organism}/seq/gsms/{gsm}.tsv",
            gsm=gsms,
            organism=wildcards.organism)
    return gsm_files

rule gsm_to_gse:
    input:
        gsm_files=postquant_gsms_for_gse,
        gse_df="out/{organism}/seq/data/filtering/postquant/gsm_gse.tsv",
        gene_mapping="input/{organism}/ensamble_genesymbol_entrez.tsv"
    output:
        gse="out/{organism}/seq/gses/{gse}.tsv"
    log: "logs/{organism}/gsm_to_gse/{gse}.log"
    message: "Aggregating {wildcards.gse} ({wildcards.organism})..."
    shadow: "shallow"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/gsm_to_gse.R {output.gse} out/{wildcards.organism}/seq/gsms {input.gene_mapping} "
        " {input.gse_df}"
        " > {log} 2>&1"

def get_postquant_passing_gse(wildcards):
    print("Getting list of GSE to be aggregated")
    filtered_gse_list = checkpoints.postquant_filter.get(**wildcards).output.passing_gse_list
    gses = [line.rstrip('\n') for line in open(filtered_gse_list)]
    gse_files = \
        expand("out/{organism}/seq/gses/{gse}.tsv",
            gse=gses,
            organism=wildcards.organism)
    return gse_files

rule gse_checkout:
    input: get_postquant_passing_gse
    output: "out/{organism}/seq/gse_checkout_flag"
    shell: "touch {output}"

#############################################CHIP_ROOT##################################################################
rule extract_exp_mat:
    '''
    Extracts expression matrices from files and writes QC report for downstream filtering.
    gq_df: TAG	N_GSM	GPL	HAS_EXP_MAT	LOGAV	LINMAX	LOGMAX	N_GENES	HAS_NEGATIVE_VALUES	PROCESSED
    Writes empty exp table if error produced.
    '''
    resources:
        time=40
    input:
        sm="out/{organism}/chip/series_matrices/{tag}_series_matrix.txt.gz",
        gpl_dir="input/{organism}/chip/platform_annotation"
    output:
        exp_table=protected("out/{organism}/chip/exp_mat/{tag}.tsv"),
        qc_report=protected("out/{organism}/chip/exp_mat/{tag}_qc.tsv")
    log: "logs/{organism}/chip/extract_exp_mat/{tag}.log"
    message: "Extracting expression table from {wildcards.tag}_series_matrix.txt.gz ({wildcards.organism})..."
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/get_exp_table.R {input.sm} {input.gpl_dir} {output.exp_table} {output.qc_report}"
        " > {log} 2>&1"

def get_filtered_sm_qc_files(wildcards):
    gse_df_file = checkpoints.extract_sm_metadata.get(**wildcards).output.gse_df
    gse_df = pd.read_csv(gse_df_file, sep="\t")
    gse_df = gse_df[(gse_df['NUMBER_GSM']>=config['min_gsm']) & (gse_df['NUMBER_GSM']<=config['max_gsm'])]
    # getting list of supported GPL
    gpl_dir = f"input/{wildcards.organism}/chip/platform_annotation"
    gpl_list = [f for f in listdir(gpl_dir) if isfile(join(gpl_dir, f))]
    gpl_list = [re.search('(.+?).3col.gz', gpl).group(1) for gpl in gpl_list]
    gse_df = gse_df[gse_df["GPL"].isin(gpl_list)]

    gse_df = gse_df[gse_df["SUPER_SERIES_OF"].isna()]

    filtered_sm_tags=gse_df["GSE"].to_list()

    filtered_sm_qc_files = \
        expand('out/{organism}/chip/exp_mat/{tag}_qc.tsv',
        organism=wildcards.organism,
        tag=filtered_sm_tags)

    return filtered_sm_qc_files

checkpoint get_exp_mat_qc_df:
    input:
        get_filtered_sm_qc_files
    output:
        "out/{organism}/{platform}/exp_qc_df.tsv"
    message: "Getting expression matrices QC df ({wildcards.organism})..."
    log: "logs/{organism}/{platform}/get_exp_mat_qc_df.log"
    run:
        qc_df_list=[pd.read_csv(qc_file,sep='\t') for qc_file in input]
        qc_df=pd.concat(qc_df_list)
        qc_df.to_csv(str(output),sep='\t',index=False)
        print("Done.")

#############################################INPUT_FUNCTION_QC##########################################################
def get_filtered_exp_mat_files(wildcards, min_gsm=int, max_gsm=int, min_genes=int,
                                       logav_min=config["chip_logav_min"], logav_max=config["chip_logav_max"],
                                       logmax_max=config["chip_logmax_max"], linmax_max=config["chip_linmax_max"],
                                       allow_negative_val=False):
    """
    QC logic and rooting for choosing files for downstream analysis.
    allow_negative_val-flag tells if negative values allowed in exp mat, which is the case for some chips.
    """
    if wildcards.platform=="chip":
        sm_qc_df_file = str(checkpoints.get_exp_mat_qc_df.get(**wildcards).output)
        sm_qc_df = pd.read_csv(sm_qc_df_file, sep="\t")
        sm_qc_df = sm_qc_df[sm_qc_df['PROCESSED']==True]
        sm_qc_df = sm_qc_df[(sm_qc_df['N_GSM']>=min_gsm) & (sm_qc_df['N_GSM']<=max_gsm)]
        sm_qc_df = sm_qc_df[(sm_qc_df['LOGAV']>=logav_min) & (sm_qc_df['LOGAV']<=logav_max)]
        sm_qc_df = sm_qc_df[sm_qc_df['LINMAX']<=linmax_max]
        sm_qc_df = sm_qc_df[sm_qc_df['LOGMAX']<=logmax_max]
        sm_qc_df = sm_qc_df[(sm_qc_df['HAS_NEGATIVE_VALUES']==False) | allow_negative_val]
        sm_qc_df = sm_qc_df[sm_qc_df['N_GENES']>=min_genes]
        exp_mat_tags = sm_qc_df["TAG"].tolist()

    elif wildcards.platform=="seq":
        gsm_gse_df_file = str(checkpoints.postquant_filter.get(**wildcards).output.gsm_gse_df)
        gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
        gse_df =  gse_df.groupby('GSE')['GSE'].transform('count')
        exp_mat_tags = gse_df[gse_df["freq"]>=min_gsm | gse_df['freq']<=max_gsm].tolist()
    return exp_mat_tags

def get_filtered_tags(wildcards, min_gsm=int, max_gsm=int, min_genes=int,
                               logav_min=config["chip_logav_min"], logav_max=config["chip_logav_max"],
                               logmax_max=config["chip_logmax_max"], linmax_max=config["chip_linmax_max"],
                               allow_negative_val=False):
    """
    QC logic and rooting for choosing files for downstream analysis.
    allow_negative_val-flag tells if negative values allowed in exp mat, which is the case for some chips.
    """
    if wildcards.platform=="chip":
        sm_qc_df_file = f"out/{wildcards.organism}/chip/exp_qc_df.tsv"
        sm_qc_df = pd.read_csv(sm_qc_df_file, sep="\t")
        sm_qc_df = sm_qc_df[sm_qc_df['PROCESSED']==True]
        sm_qc_df = sm_qc_df[(sm_qc_df['N_GSM']>=min_gsm) & (sm_qc_df['N_GSM']<=max_gsm)]
        sm_qc_df = sm_qc_df[(sm_qc_df['LOGAV']>=logav_min) & (sm_qc_df['LOGAV']<=logav_max)]
        sm_qc_df = sm_qc_df[sm_qc_df['LINMAX']<=linmax_max]
        sm_qc_df = sm_qc_df[sm_qc_df['LOGMAX']<=logmax_max]
        sm_qc_df = sm_qc_df[(sm_qc_df['HAS_NEGATIVE_VALUES']==False) | allow_negative_val]
        sm_qc_df = sm_qc_df[sm_qc_df['N_GENES']>=min_genes]
        exp_mat_tags = sm_qc_df["TAG"].tolist()

    elif wildcards.platform=="seq":
        gsm_gse_df_file = str(f"out/{wildcards.organism}/seq/postquant_filter/gsm_gse.tsv")
        gse_df = pd.read_csv(gsm_gse_df_file, sep="\t")
        gse_df =  gse_df.groupby('GSE')['GSE'].transform('count')
        exp_mat_tags = gse_df[gse_df["freq"]>=min_gsm | gse_df['freq']<=max_gsm].tolist()

    return exp_mat_tags#[:len(exp_mat_tags)//40]


#############################################PCA########################################################################
rule pca:
    '''
    Takes exp_mat as an input and outputs first 
    1) 10 (or less) PC components loadings scores
    2) Stats- explained var, etc. 
    '''
    resources:
        time = 20
    input:
        "out/{organism}/{platform}/exp_mat/{tag}.tsv"
    output:
        protected("out/{organism}/{platform}/pca/{max_genes}_{scale}/{tag}.rds"),
    message:
        "Performing PCA on {wildcards.tag} \n"
        " Organism: {wildcards.organism} \n"
        " Platform: {wildcards.platform} \n"
        " Number of genes considered: {wildcards.max_genes} \n"
        " Scale: {wildcards.scale}"
    log: "logs/{organism}/{platform}/pca/{max_genes}_{scale}/{tag}.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/pca.R {input} {output} {wildcards.max_genes} {wildcards.scale}"
        " > {log} 2>&1"

rule get_pc_list_adapter:
    '''
    Apparently it is impossible to pass 4k arguments with bash. This is workaround with concatenation of arguments to a 
    single list. Note: it is not a problem with large number of inputs if they passed to 'run' section.
    '''
    input:
        lambda wildcards:
            expand(rules.pca.output,
                organism=wildcards.organism,
                platform=wildcards.platform,
                max_genes=wildcards.max_genes,
                scale=wildcards.scale,
                tag=get_filtered_exp_mat_files(wildcards, int(config["pca_min_gsm"]),
                int(config["pca_max_gsm"]), int(wildcards.max_genes)))
    output:
        "out/{organism}/{platform}/pca/{max_genes}_{scale}.list"
    run:
        with open(str(output), 'w') as f:
            for file in input:
                f.write("%s\n" % str(file))

rule get_pc_list:
    '''
    Aggregates pca components to single list of ranked lists. Filters components by QC.
    QC params:
    1) Explained valiance threshold. 2% for starters.
    2) No more then 10pc per df.
    If results look good conditions might be relaxed.
    '''
    resources:
        mem_ram=20
    input:
        rules.get_pc_list_adapter.output
    message:
        "Getting PC list {wildcards.organism} {wildcards.platform} \n"
        " Number of genes considered: {wildcards.max_genes} \n"
        " Scale of original dataset: {wildcards.scale} \n"
        " Explained variance % threshold: {wildcards.var_threshold} \n"
        " Max PC components for 1 dataset: {wildcards.max_comp} \n"
    log: "logs/{organism}/{platform}/get_pc_list/{max_genes}_{scale}_{max_comp}_{var_threshold}.log"
    output:
        "out/{organism}/{platform}/pca/{max_genes}_{scale}_{max_comp}_{var_threshold}_PCList.rds"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/get_pc_list.R {output} {wildcards.max_comp} {wildcards.var_threshold} {input}"
        " > {log} 2>&1"

# TODO: remove first two lines in genesets inside the script. Artifact from GQ
rule fgsea_genesets:
    '''
    Performs fgsea against list of PC components. Outputs ranked list of results with NES.
    '''
    resources:
        mem_ram=20,
        time=360,
        cores=8
    input:
        pc_list=rules.get_pc_list.output,
        geneset="input/{organism}/genesets/{geneset_name}",
    output:
        "out/{organism}/{platform}/pca_fgsea/"
        "{max_genes}_{scale}_{max_comp}_{var_threshold}/"
        "raw/{geneset_name}.tsv"
    message:
        "Performing GSEA {wildcards.organism} {wildcards.platform} \n"
        " Geneset: {wildcards.geneset_name} \n"
        " Number of genes considered: {wildcards.max_genes} \n"
        " Scale of original dataset: {wildcards.scale} \n"
        " Explained variance threshold %: {wildcards.var_threshold} \n"
        " Max PC components for 1 dataset: {wildcards.max_comp}"
    log:
        "logs/{organism}/{platform}/fgsea_genesets/"
        "{max_genes}_{scale}_{max_comp}_{var_threshold}/"
        "{geneset_name}.log"
    conda: "envs/fgsea.yaml"
    shell:
        "Rscript scripts/R/fgsea_geneset.R {input.pc_list} {input.geneset} {output}"
        " > {log} 2>&1"

rule prepare_pca_fgsea_result:
    input:
        gsea_results=rules.fgsea_genesets.output,
        gse_df="out/{organism}/{platform}/sm_metadata/gse.tsv"
    output:
        "out/{organism}/{platform}/pca_fgsea/{max_genes}_{scale}_{max_comp}_{var_threshold}/"
        "prepared/{geneset_name}.tsv"
    message:
        "Preparing results for PCA query. {wildcards.organism} {wildcards.platform} \n"
        " Geneset: {wildcards.geneset_name} \n"
        " Number of genes considered: {wildcards.max_genes} \n"
        " Scale of original dataset: {wildcards.scale} \n"
        " Explained variance threshold %: {wildcards.var_threshold} \n"
        " Max PC components for 1 dataset: {wildcards.max_comp}"
    log:
        "logs/{organism}/{platform}/prepare_pca_fgsea_result/"
        "{max_genes}_{scale}_{max_comp}_{var_threshold}/"
        "{geneset_name}.log"
    conda: "envs/r_scripts.yaml"
    shell:
        "Rscript scripts/R/pca_prepare_results.R {input.gsea_results} {input.gse_df} {output}"
        " > {log} 2>&1"


#########################################WGCNA##########################################################################
# # TODO: CHANGE GSES in seq to exp_mat so that it is assessable from downstream rules with similar path pattern
# rule wgcna:
#     input:
#         "out/{organism}/{platform}/exp_mat/{tag}.tsv"
#     output:
#         modules="out/{organism}/{platform}/wgcna/{n_genes}_{cor_type}_{network}_{scale}/{tag}/modules.tsv",
#         eigengenes="out/{organism}/{platform}/wgcna/{n_genes}_{cor_type}_{network}_{scale}/{tag}/eigengenes.tsv",
#         stats="out/{organism}/{platform}/wgcna/{n_genes}_{cor_type}_{network}_{scale}/{tag}/stats.txt"
#     log: "logs/{organism}/{platform}/wgcna/{n_genes}_{cor_type}_{network}_{scale}/{tag}.log"
#     message: "Performing WGCNA on {wildcards.tag} ({wildcards.organism}, {wildcards.platform})..."
#     conda: "envs/r_scripts.yaml"
#     shell:
#         "Rscript scripts/R/wgcna.R {input} {output.modules} {output.eigengenes} {output.stats} "
#         " {wildcards.n_genes} {wildcards.cor_type} {wildcards.network} {wildcards.scale}"
#         " > {log} 2>&1"
# #
# checkpoint get_wgcna_stats:
#     input:
#         lambda wildcards:
#             expand("out/{organism}/{platform}/wgcna/{n_genes}_{cor_type}_{network}_{scale}/{tag}/stats.txt",
#                 organism=wildcards.organism,
#                 platform=wildcards.platform,
#                 n_genes=wildcards.n_genes,
#                 cor_type=wildcards.cor_type,
#                 network=wildcards.network,
#                 scale=wildcards.scale,
#                 tag=get_filtered_exp_mat_files(wildcards, int(config["min_gsm_wgcna"]),
#                 int(config["max_gsm_wgcna"]), int(wildcards.n_genes)))
#     output:
#         "out/{organism}/{platform}/wgcna/{n_genes}_{cor_type}_{network}_{scale}_stats.tsv"
#     conda: "envs/r_scripts.yaml"
#     log: "logs/{organism}/{platform}/get_wgcna_stats/{n_genes}_{cor_type}_{network}_{scale}.log"
#     message:
#         "Getting WGCNA stats. \n"
#         " Organism: {wildcards.organism} \n"
#         " Platform: {wildcards.platform} \n"
#         " Number of genes considered: {wildcards.n_genes} \n"
#         " Correlation metric: {wildcards.cor_type} \n"
#         " Network type: {wildcards.network} \n"
#         " Scale: {wildcards.scale} \n"
#     run:
#         qc_df_list=[pd.read_csv(qc_file,sep='\t') for qc_file in input]
#         qc_df=pd.concat(qc_df_list)
#         qc_df.to_csv(str(output),sep='\t',index=False)
#         print("Done.")
