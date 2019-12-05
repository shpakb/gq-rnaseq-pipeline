#!/bin/bash
#####################################################################################1
# Takes geo data sets search result in .txt format
# Outputs dir with series matrices
gdsResult=$1
outDir=$2

echo "GEO search results file: $gdsResult"
echo "Out dir: $outDir"

# extract URL prefixes from the file
URLS=$(grep "FTP download" "$gdsResult" | grep GSE | perl -ne 'm/GEO( | \(.*\) )(ftp\:\/\/ftp.*\/)/; print "$2\n"' | awk -F "," '{print $1}')

for i in $URLS; do
  {
    echo "Downloading from address $i"
    wget -nc -nv -P "$outDir" "$i"/matrix/*_series_matrix.txt.gz
  } ||
  {
    echo "Trying again..."
    wget -nc -nv -P "$outDir" "$i"/matrix/*_series_matrix.txt.gz
  } ||
  {
    echo "Failed"
  }
done

echo "All downloads are now complete!"
