#!/usr/bin/env python3

#Takes an unsorted, unindexed bam file
#Uses samtools to sort, and index it
#Uses sinto to fragment it
#Uses samtools to sort, zip, and index the fragments
import argparse
from argparse import ArgumentParser
from subprocess import check_call
import os

SORTED_FILE_NAME = "sorted_bam_file.bam"
BED_FILE_NAME = "frags.bed"
SORTED_BED_FILE = "frags.sort.bed"
ZIPPED_BED_FILE = "frags.sort.bed.gz"

#Steps in pipeline that don't require runtime input
steps = ["index_bam", "generate_fragments", "sort_bed", "zip_bed", "index_bed"]

#Mapping from step of pipeline to command for check_call
commands = {"index_bam": ["samtools", "index", SORTED_FILE_NAME],
"generate_fragments": ["sinto", "fragments", "-b", SORTED_FILE_NAME, "-f", BED_FILE_NAME],
"sort_bed": ["sort", "-k1,1", "-k2,2n", BED_FILE_NAME, ">", SORTED_BED_FILE],
"zip_bed": ["bgzip", SORTED_BED_FILE],
"index_bed": ["tabix", "-p", "bed", ZIPPED_BED_FILE]}

def sort_bam(bam_file: str):
    print("Bam sorting")
    command = ["samtools", "sort", bam_file, "-o", SORTED_FILE_NAME]
    check_call(command)
    return

def main(bam_file: str):
    #This step requires runtime input, so it's separate
    print("Program executed")
    sort_bam(bam_file)
    print("Bam sorted")
    for step in steps:
        print(step)

        if step == "sort_bed":#This is super hacky, but check_call chokes when we redirect input, but system works
            this_command = ""
            for command in commands[step]:
                this_command += command
                this_command += " "
            os.system(this_command)

        else:
           check_call(commands[step])
    #@TODO:Cleanup sorted_bam_file, frags.bed, etc?

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bam_file', type=str)
    args = p.parse_args()

    main(args.bam_file)
