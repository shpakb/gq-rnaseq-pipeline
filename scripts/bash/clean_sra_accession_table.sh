#!/bin/bash

inFile=$1
outFile=$2

echo "Input file: $inFile"
echo "Output file: $outFile"

echo "SRR GSM SPOTS" > "$outFile"

awk '{print $1,$10,$15}' "$inFile" | grep SRR | grep GSM >> "$outFile"

echo "Done."
