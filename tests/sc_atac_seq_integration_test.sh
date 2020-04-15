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

   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/BUKMAP_20190529I_1000000__R1.fastq
   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/BUKMAP_20190529I_1000000__R2.fastq
   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/BUKMAP_20190529I_1000000__R3.fastq

   cwltool --debug --timestamps --target snap_file $CWLTOOL_TMPDIR_PREFIX $CWLTOOL_TMP_OUTDIR_PREFIX  ../create_snap_and_analyze.cwl create_snap_and_analyze.json
  


