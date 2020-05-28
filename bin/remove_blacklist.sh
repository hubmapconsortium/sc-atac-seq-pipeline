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

#Set default blacklist file
encode_blacklist=/opt/supplementary-data/hg38.blacklist.bed
#Overwrite default if a blacklist file is provided
if [[ "$*" == *--bed-file* ]]
 then
   encode_blacklist=${4}
fi

#If the user provided a blacklist file, or is using the default alignment_index, remove blacklist
if [[ "$*" == *--bed-file* ]] || [[ "$*" != *--alignment-index* ]]
then
  echo "Removing blacklisted regions"
  # Get the PG tags from the header
  # SnapTools snapPre needs thisgrch38_index.tar.gz
  # and bedtools intersect strips it from the header
  samtools view -H ${2} | grep @PG > pg.txt

  bedtools intersect -v -abam ${2} -b $encode_blacklist > unsorted_rmsk.bam
  # The resultant BAM file will be unsorted.
  # Sort the resulting bam file by name for the next step in the pipeline
  samtools sort -n -o no_pg_rmsk.bam -T tmp unsorted_rmsk.bam

  # Extract the new header from the BAM file
  samtools view -H no_pg_rmsk.bam > no_pg_rmsk_header.txt
  # Concatenate the PG tag to the end of the header
  cat pg.txt >> no_pg_rmsk_header.txt
  # Use the new header in the output BAM file
  samtools reheader no_pg_rmsk_header.txt no_pg_rmsk.bam > rmsk.bam
#If the user provided an alignment_index but not a blacklist file, don't remove anything
else
  cp ${2} rmsk.bam
fi
