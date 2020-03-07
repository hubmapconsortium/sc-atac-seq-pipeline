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

# https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
filename1nopath=$(basename -- "${1}")
extension1="${filename1nopath#*.}"
filename1="${filename1nopath%+(.fastq.gz|.fastq)}"

filename2nopath=$(basename -- "${2}")
extension2="${filename2nopath#*.}"
filename2="${filename2nopath%%+(.fastq.gz|.fastq)}"

cat  "${1}" "${2}" > "${filename1}"_"${filename2}"."${extension1}"

