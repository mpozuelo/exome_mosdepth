#!/usr/bin/env bash

folder=$1
output=$2

for f in $(find $folder -name "*.txt")
do
name=$(echo $f | sed 's/\.\/\([A-Za-z0-9_]*\)\/[A-Za-z0-9_\.]*/\1/g')
stats=$(echo $(grep -E "^BAIT" -A 1 $f | grep -v BAIT))
echo $name $stats >> $output
done

header=$(echo $(grep -E ^BAIT $f))
complete_header=$(echo "sampleID" $header)
sed -i "1s/^/$complete_header\n/g" $output
