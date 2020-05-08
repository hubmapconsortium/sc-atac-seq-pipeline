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

   # Download these files using wget because when https links were provided in the JSON input file
   # the download threw a memory error.
#   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/grch38.fasta
#   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/grch38_index.tar.gz
#   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/hg38.blacklist.bed
#   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/gencode.v32.annotation.bed
#   wget https://storage.googleapis.com/sc-atac-seq-pipeline-testing/hg38.promoters.bed

   cwltool --debug --timestamps $CWLTOOL_TMPDIR_PREFIX $CWLTOOL_TMP_OUTDIR_PREFIX  ../create_snap_and_analyze.cwl create_snap_and_analyze_local.json
  


