#!/bin/bash
  
# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit 
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u
# echo each line of the script to stdout so we can see what is happening
set -o xtrace
#to turn off echo do 'set +o xtrace'

input_sam=$1
# Get the cell id from the fastq comment which is the part of the comment before the first colon
# Add the CB (cell identifier tag) to the end of the read using the cell id
# https://stackoverflow.com/questions/19323529/how-can-i-get-sed-to-only-match-to-the-first-occurrence-of-a-character
# https://stackoverflow.com/questions/9053100/sed-regex-and-substring-negation
sed -E '/^@/!s/^([^:]*).*/&\tCB:Z:\1/g' "${input_sam}" > CB.sam
# Convert the SAM file to a BAM file
samtools sort -O sam -T sample.sort -o CB_sorted.sam CB.sam
samtools view -S -b CB_sorted.sam > CB_sorted.bam

