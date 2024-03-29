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

if [[ "$*" == *--bed-file* ]]
then
  echo "Removing blacklisted regions"
  # Get the PG tags from the header
  # TODO is it necessary to do this anymore? SnapTools snapPre needed this
  # and bedtools intersect strips it from the header
  samtools view -H ${2} | grep @PG > pg.txt
  
  bedtools intersect -v -abam ${2} -b ${4} > unsorted_rmsk.bam
  # The resultant BAM file will be unsorted.
  # Sort the resulting bam file by name for the next step in the pipeline
  samtools sort -n -o no_pg_rmsk.bam -T tmp unsorted_rmsk.bam
  
  # Extract the new header from the BAM file
  samtools view -H no_pg_rmsk.bam > no_pg_rmsk_header.txt
  # Concatenate the PG tag to the end of the header
  cat pg.txt >> no_pg_rmsk_header.txt
  # Use the new header in the output BAM file
  samtools reheader no_pg_rmsk_header.txt no_pg_rmsk.bam > rmsk.bam
else
  cp ${2} rmsk.bam
  echo "Skipping blacklist removal; no genome BED file with the blacklisted regions to be removed provided"
fi


