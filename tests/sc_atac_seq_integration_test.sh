#! /bin/bash

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

   pwd
   #docker version
   #cwltool --version

   #ls -al
   #df -h
   #sudo du -hsx ./* | sort -n | head -100
   #export $TMPDIR=/mnt/tmp
   #set


   cwltool --debug --timestamps --target analysis_TXT_files $CWLTOOL_TMPDIR_PREFIX $CWLTOOL_TMP_OUTDIR_PREFIX  ../steps/create_snap_and_analyze.cwl create_snap_and_analyze.json



