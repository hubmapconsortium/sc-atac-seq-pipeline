#!/usr/bin/env pypy3
import shlex
from argparse import ArgumentParser
from multiprocessing import Pool
from pathlib import Path
from subprocess import PIPE, Popen, check_call
from typing import List

ALIGNED_BAM_FILENAME = "alignment.bam"

samtools_sort_command_template = [
    "samtools",
    "sort",
    "--threads",
    "{processes}",
    "-o",
    ALIGNED_BAM_FILENAME,
    "{merged_bam}",
]
samtools_index_command_template = [
    "samtools",
    "index",
    "-@",
    "{processes}",
    ALIGNED_BAM_FILENAME,
]


def adjust_sam_line(line: bytes) -> bytes:
    barcode = line.split(b":")[0]
    return line.strip() + b"\tCB:Z:" + barcode + b"\n"


def print_command(command: List[str]):
    print("Running", " ".join(shlex.quote(c) for c in command))


def main(processes: int, merged_bam_file: Path):

    sort_and_index_commands = [
        [piece.format(processes=processes, merged_bam=merged_bam_file) for piece in samtools_sort_command_template],
        [piece.format(processes=processes) for piece in samtools_index_command_template],
    ]
    for command in sort_and_index_commands:
        print_command(command)
        check_call(command)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("-p", "--processes", type=int, default=1)
    p.add_argument("merged_bam_file", type=Path)
    args = p.parse_args()

    main(args.processes, args.merged_bam_file)
