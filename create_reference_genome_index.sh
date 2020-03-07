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

# If there is a reference index tar gz supplied
if [[ "$*" == *--reference-genome-index* ]]
then
  echo "Untarring reference genome index tar gz"
  while test $# -gt 0
  do
    case "$1" in
        --reference-genome-index)
            echo "Untarring reference genome index tar gz file $1"
            shift
            tar -xvzf $1
            exit 0;
            ;;
        *) echo "Error did not find reference genome argument"
            ;;
    esac
    shift
  done
else
  echo "Indexing reference genome"
  snaptools index-genome $*
fi
