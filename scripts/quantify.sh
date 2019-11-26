#!/bin/bash

SRR=$1
FASTQ_DIR=$2
REFSEQ=$3
OUTPUT_DIR=$4

SRR_PREF=$FASTQ_DIR/$SRR

echo "SRR: $SRR"
echo "Fastq dir: $FASTQ_DIR"
echo "REFSEQ: $REFSEQ"
echo "Output dir: $OUTPUT_DIR"
echo "Working directory: $(pwd)"

N=$(find $FASTQ_DIR -wholename "$SRR_PREF*fastq" | wc -l )

if [[ $N == 2 ]]
then
  echo "Two fastq files found; processing sample $SRR as a paired-ended experiment."
  echo "$REFSEQ -o $SRR $SRR_PREF*.fastq"
  kallisto quant -i $REFSEQ -o $OUTPUT_DIR $SRR_PREF*.fastq
elif [[ $N == 3 ]]
then
  echo "Three fastq files found; removing single-end reads and processing sample $SRR as a paired-ended experiment."
  echo "$REFSEQ -o $SRR $SRR_PREF*.fastq"
  kallisto quant -i $REFSEQ -o $OUTPUT_DIR $SRR_PREF*.fastq
elif [[ $N == 1 ]]
then
  echo "One fastq file found; processing sample $SRR as a single-ended experiment."
  echo "$REFSEQ -o $SRR $SRR_PREF*.fastq"
  kallisto quant --single -l 200 -s 50 -i $REFSEQ -o $OUTPUT_DIR $SRR_PREF*.fastq
else 
  echo "ERROR: Wrong number of input arguments!"
  exit 1
fi

echo "Quantification complete."