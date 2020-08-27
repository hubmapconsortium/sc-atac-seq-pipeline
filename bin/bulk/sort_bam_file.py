#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import check_call

SORT_COMMAND_TEMPLATE = [
    "samtools",
    "sort",
    "{bam_file}",
    "-o",
    "{sorted_bam_name}",
    "--threads",
    "{threads}",
]

def main(bam_file: Path, sorted_bam_name: Path, threads: int):
    command = [
        piece.format(
            bam_file=bam_file,
            sorted_bam_name=sorted_bam_name,
            threads=threads,
        )
        for piece in SORT_COMMAND_TEMPLATE
    ]
    check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bam_file', type=Path)
    p.add_argument('sorted_bam_name', type=Path)
    p.add_argument('--threads', type=int, default=0)
    args = p.parse_args()

    main(args.bam_file, args.sorted_bam_name, args.threads)
