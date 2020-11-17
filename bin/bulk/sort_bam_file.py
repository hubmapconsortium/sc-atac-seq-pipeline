#!/usr/bin/env python3
"""
Takes an unsorted, unindexed bam file
Uses samtools to sort, and index it
Uses sinto to fragment it
Uses samtools to sort, zip, and index the fragments
"""
import argparse
import os
from argparse import ArgumentParser
from subprocess import check_call


def main(bam_file: str, sorted_bam_name: str):
    command = ["samtools", "sort", bam_file, "-o", sorted_bam_name]
    check_call(command)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("bam_file", type=str)
    p.add_argument("sorted_bam_name", type=str)
    args = p.parse_args()

    main(args.bam_file, args.sorted_bam_name)
