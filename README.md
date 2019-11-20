# GQ rnaseq quantification pipeline

Pipeline quantifies SRR files and aggregates them first to GSM and then to GSM matrices.


### 1) Run docker container:
- Better run in tmux as interactive session:
```bash
bsub -Is  -q docker-interactive  -a 'docker(biolabs/snakemake:5.6.0_conda4.7.12)' /bin/bash
```
### 2) Git clone snakemake repo to home directory:
```bash
git clone https://github.com/shpakb/gq-rnaseq-pipeline.git
```

### 3) Create storage pipeline_wdir at storage location 
```bash 
mkdir /gscmnt/gc2676/martyomov_lab/shpakb/pipeline_wdir && \
    cd /gscmnt/gc2676/martyomov_lab/shpakb/pipeline_wdir
```
### 4) Create symlinks and copying necessary stuff from pipeline repo to pipeline_wdir :

- inside pipline_wdir
```bash
ln -s /gscmnt/gc2676/martyomov_lab/shpakb/Assemblies/rnor_v6/ index && \
    ln -s ~/gq-rnaseq-pipeline/envs/ envs && \
    ln -s ~/gq-rnaseq-pipeline/scripts/ scripts && \
    cp ~/gq-rnaseq-pipeline/config.yaml . && \
    cp ~/gq-rnaseq-pipeline/lsf_jobscript.sh . && \
    cp ~/gq-rnaseq-pipeline/Snakefile .
```

### 5) Put srr.list file with that has to be quantified to pipline_wdir

### 6) Create and activate conda env
```bash 
conda env create --file ./envs/quantify.yaml --name snakemake && \
    source activate snakemake
```

### 6) Run pipeline:
```bash
snakemake -pr --use-conda --profile lsf --jobs 50 --restart-times 3 \
    --directory $(pwd) \
    --jobscript $(pwd)/lsf_jobscript.sh
    --resources load=100
```

- resources parameter specifies amount of resources pipline can use. In this particular case load 100 and 
each downloading job uses 50 "points" of load. So ther can't be more than two downloading jobs simultaneously. 

