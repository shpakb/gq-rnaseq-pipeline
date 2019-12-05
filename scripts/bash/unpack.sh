#!/bin/bash

OUTDIR=/gscmnt/gc2676/martyomov_lab/shpakb/gq-rnaseq-pipeline/out/rn/seq/kallisto
INDIR=/gscmnt/gc2676/martyomov_lab/shpakb/rna_seq/all_rn
SRRS=$( ls ${INDIR}/abund | sed -e 's/\.abundance.tsv$//')


for SRR in SRRS; do
    echo $SRR
    {
        SRRDIR=${OUTDIR}/{SRR}
        mkdir $SRRDIR
        cp ${INDIR}/abund/${SRR}.abundance.tsv ${SRRDIR}/abundance.tsv
        cp ${INDIR}/h5/${SRR}.abundance.h5 ${SRRDIR}/abundance.h5
        cp ${INDIR}/json/${SRR}.run_info.json ${SRRDIR}/run_info.json
    } ||
    {
        echo "Missing one of the files"
        rm -rf $SRRDIR
    }
done