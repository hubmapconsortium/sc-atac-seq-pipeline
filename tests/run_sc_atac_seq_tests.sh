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

#https://unix.stackexchange.com/questions/129391/passing-named-arguments-to-shell-scripts
while getopts ":t:p:o:" options; do
    case "$options" in
      t )
        echo "found t"
        CWLTOOL_TMPDIR="$OPTARG"
        ;;
      p )
        echo "found p"
        CWLTOOL_TMPDIR_PREFIX="--tmpdir-prefix $OPTARG"
        ;;
      o )
         echo "found o"
         CWLTOOL_TMP_OUTDIR_PREFIX="--tmp-outdir-prefix $OPTARG"
         ;;
      \? ) 
         echo "Invalid option -$OPTARG" >&2
         exit 1 
         ;;
      : )
         echo "Invalid option: $OPTARG requires an argument" 1>&2
         ;;
    esac
done
shift $((OPTIND -1))

## cwltool writes files downloaded using http or https to the directory specified
## by the TMPDIR env var 
#https://unix.stackexchange.com/questions/122845/using-a-b-for-variable-assignment-in-scripts
TMPDIR="${TMPDIR:-$CWLTOOL_TMPDIR}"
CWLTOOL_TMPDIR_PREFIX="${CWLTOOL_TMPDIR_PREFIX:-}"
CWLTOOL_TMP_OUTDIR_PREFIX="${CWLTOOL_TMP_OUTDIR_PREFIX:-}"

echo "$TMPDIR"
echo "$CWLTOOL_TMPDIR_PREFIX"
echo "$CWLTOOL_TMP_OUTDIR_PREFIX"

pwd
# cd to directory where test scripts are located
cd $(dirname $0)
pwd
# https://stackoverflow.com/questions/8352851/how-to-call-one-shell-script-from-another-shell-script
# using source executes the script in the first script's process, and pulls in
# variables and functions from the other script so they are usable from the calling script.
source ./sc_atac_seq_integration_test.sh
