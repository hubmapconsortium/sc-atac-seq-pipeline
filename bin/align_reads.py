#!/usr/bin/env pypy3
import shlex
from argparse import ArgumentParser
from multiprocessing import Pool
from pathlib import Path
from subprocess import PIPE, Popen, check_call
from typing import List

align_command_template = [
    "hisat2",
    "-x",
    "/opt/supplementary-data/hisat2-index/hg38/genome",
    "-p",
    "{processes}",
    "-1",
    "{fastq_1}",
    "-2",
    "{fastq_2}",
]
samtools_view_command = [
    "samtools",
    "view",
    "-h",
    "-b",
    "-o",
    "unsorted.bam",
    "-",
]
samtools_sort_command_template = [
    "samtools",
    "sort",
    "--threads",
    "{processes}",
    "-o",
    "alignment.bam",
    "unsorted.bam",
]


def adjust_sam_line(line: bytes) -> bytes:
    barcode = line.split(b":")[0]
    return line.strip() + b"\tCB:Z:" + barcode + b"\n"


def print_command(command: List[str]):
    print("Running", " ".join(shlex.quote(c) for c in command))


def main(processes: int, fastq_1: Path, fastq_2: Path):
    align_command = [
        piece.format(processes=processes, fastq_1=fastq_1, fastq_2=fastq_2)
        for piece in align_command_template
    ]
    print_command(align_command)
    align_proc = Popen(align_command, stdout=PIPE)

    print_command(samtools_view_command)
    samtools_view_proc = Popen(samtools_view_command, stdin=PIPE)
    for line in align_proc.stdout:
        if line.startswith(b"@"):
            samtools_view_proc.stdin.write(line)
        else:
            samtools_view_proc.stdin.write(adjust_sam_line(line))
            break
    # rest of lines need to be adjusted now
    with Pool(processes // 2) as pool:
        for adj_line in pool.imap_unordered(adjust_sam_line, align_proc.stdout):
            samtools_view_proc.stdin.write(adj_line)

    samtools_view_proc.stdin.close()
    samtools_view_proc.wait()
    align_proc.wait()

    samtools_sort_command = [
        piece.format(processes=processes) for piece in samtools_sort_command_template
    ]
    print_command(samtools_sort_command)
    check_call(samtools_sort_command)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("-p", "--processes", type=int, default=1)
    p.add_argument("fastq_1", type=Path)
    p.add_argument("fastq_2", type=Path)
    args = p.parse_args()

    main(args.processes, args.fastq_1, args.fastq_2)
