#!/bin/bash

inFile=$1
outFile=$2

echo "Input file: $inFile"
echo "Output file: $outFile"

RES=$(cut -f 1,10 "$inFile" | perl -ne '@t=split/\t/; $t[1]=~m/(GSM\d+)/; print "$t[0]\t$1\n"')

RES=$(grep ^SRR "$RES" | grep GSM)

# Extract SRX, GSM, SPOTS columns. SPOTS is number of reads
RES=$(awk '{print $1,$10,$15}' "$RES")

$RES > outFile

echo "Done."