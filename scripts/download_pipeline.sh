#!/bin/bash
WDIR=$(pwd)

SRR=$1
OUTPUT=$2

echo "Input file: $SRR"
echo "Output file: $OUTPUT"
echo "Working directory: $WDIR"

MASK=$(echo $SRR | perl -ne 'printf "%s\n",substr $_,0,6')
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$MASK/$SRR/$SRR.sra
if [ -f "$SRR.sra" ]; then
  mv ${SRR}.sra $OUTPUT
  echo "Downloaded with wget"
else
  echo "No file in FTP directory. Using prefetch"
  prefetch $SRR --output-directory $WDIR
  mv ${SRR}/${SRR}.sra $OUTPUT
  rm -r ${SRR}
fi

echo "Done."
