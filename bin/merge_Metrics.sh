#!/usr/bin/env bash

folder=$1
output=$2
output1=$(echo $output | sed 's/_tmp.txt/.txt/g')

for f in $(find $folder -name "*.txt")
do
name=$(echo $f | sed 's/\.\/\([A-Za-z0-9_]*\)\/[A-Za-z0-9_\.]*/\1/g')
stats=$(echo $(grep -E "^BAIT" -A 1 $f | grep -v BAIT))
echo $name $stats >> $output
done

header=$(echo $(grep TOTAL_READS $f))
complete_header=$(echo "sampleID" $header)
echo $complete_header > $output1
cat $output >> $output1
