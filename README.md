# GQ rnaseq quantification pipeline

Pipeline quantifies SRR files and aggregates them first to GSM and then to GSM matrices.

### 1) Run docker container:
- Better run in tmux as interactive session:
```bash
bsub -Is  -q docker-interactive  -a 'docker(biolabs/snakemake:5.8.1_conda4.7.12)' /bin/bash
```
### 2) Git clone snakemake repo to home directory:
```bash
git clone https://github.com/shpakb/gq-rnaseq-pipeline.git
```

### 3) Create input dir with:
### - Symlinks for reference genomes.
### - Search results from GEO to parse 
### see the config file for rest 

Example how to create symlink:
```bash
ln -s /gscmnt/gc2676/martyomov_lab/shpakb/Assemblies/rnor_v6/ index 
```

### 6) Create and activate conda env(when loacal) 
```bash 
conda env create --file ./envs/quantify.yaml --name snakemake && \
    source activate snakemake
```

Test run: 
```bash
snakemake -pr --notemp --use-conda --verbose
```
Dry run(graph eval):
```bash
snakemake --dryrun
```

### 6) Run pipeline in test mode on cluster(when on cluster): 
```bash
snakemake -pr --use-conda --profile lsf --jobs 50\
--jobscript lsf_jobscript.sh\
 --resources download_res=4 writing_res=20\
 --verbose
```

Add when script is more or less stable:
add: --restart-times 3 
remove: --notemp 
--verbose
 -pr- gives command arguments overview

- resources parameter specifies amount of resources pipeline can use. In this particular case load 100 and 
each downloading job uses 50 "points" of load. So they can't be more than two downloading jobs simultaneously. 

- now pipeline outputs all the files right in to directory with scripts. No need to go through all steps with symlinks 
and copying. 


Run test snakemake:
snakemake --snakefile=test.smk