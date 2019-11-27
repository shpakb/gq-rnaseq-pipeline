#!/bin/bash

inFile=$1
outFile=$2
TEMP1=$3
TEMP2=$4

echo "Input file: $inFile"
echo "Output file: $outFile"

cut -f 1,10,14 "$inFile" | perl -ne '@t=split/\t/; $t[1]=~m/(GSM\d+)/; print "$t[0]\t$1\n"' > "$TEMP1"

grep ^SRR "$TEMP1" | grep GSM > "$TEMP2"

# Extract SRX, GSM, SPOTS columns. SPOTS is number of reads
awk '{print $1,$10,$15}' "$TEMP2" > "$outFile"

echo "Done."

