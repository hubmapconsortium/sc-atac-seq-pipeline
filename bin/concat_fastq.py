#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from pprint import pprint
from subprocess import run
from typing import Iterable, Tuple

FASTQ_EXTENSIONS = [
    "fastq",
    "fastq.gz",
]

# No point in involving the Python interpreter in this decompression;
# it'll be a lot faster just to pipe 'gunzip -c' to each file
GUNZIP_COMMAND_TEMPLATE = [
    "gunzip",
    "-c",
    "{fastq_file}",
]

FOUND_PAIR_COLOR = "\033[01;32m"
UNPAIRED_COLOR = "\033[01;31m"
NO_COLOR = "\033[00m"

MERGED_FASTQ_R1 = Path("merged_R1.fastq")
MERGED_FASTQ_R2 = Path("merged_R2.fastq")
MERGED_FASTQ_BARCODE = Path("merged_barcode.fastq")


def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = "**/*_R1*.{extension}"
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))


def find_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path, Path]]:
    """
    :param directory:
    :return: Iterable of 3-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
     [2] R3 FASTQ file
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        r2_fastq_filename = r1_fastq_file.name.replace("_R1", "_R2")
        r3_fastq_filename = r1_fastq_file.name.replace("_R1", "_R3")
        r2_fastq_file = r1_fastq_file.with_name(r2_fastq_filename)
        r3_fastq_file = r1_fastq_file.with_name(r3_fastq_filename)

        if r2_fastq_file.is_file() and r3_fastq_file.is_file():
            print(FOUND_PAIR_COLOR + "Found three FASTQ files:" + NO_COLOR)
            print(f"\t{r1_fastq_file}")
            print(f"\t{r2_fastq_file}")
            print(f"\t{r3_fastq_file}")
            yield r1_fastq_file, r2_fastq_file, r3_fastq_file
        else:
            print(UNPAIRED_COLOR + "Found lone FASTQ file:" + NO_COLOR)
            print("\t", r1_fastq_file, sep="")


def decompress_concat_fastq(input_fastq: Path, merged_fastq: Path):
    print("Concatenating", input_fastq, "to", merged_fastq)
    with open(merged_fastq, "ab") as out:
        command = [piece.format(fastq_file=input_fastq) for piece in GUNZIP_COMMAND_TEMPLATE]
        run(command, stdout=out, check=True)


def main(directories: Iterable[Path]):
    for directory in directories:
        for r1_fastq, r2_fastq, r3_fastq in find_fastq_files(directory):
            # Flipping R2 and R3 here is deliberate. For the input data, barcodes
            # are in R2, and R1 and R3 are the two reads which should be aligned
            # to the genome. Internally in this pipeline, it's more straightforward
            # to rename these: R1 -> R1, R2 -> barcode, R3 -> R2
            decompress_concat_fastq(r1_fastq, MERGED_FASTQ_R1)
            decompress_concat_fastq(r2_fastq, MERGED_FASTQ_BARCODE)
            decompress_concat_fastq(r3_fastq, MERGED_FASTQ_R2)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path, nargs="+")
    args = p.parse_args()

    main(args.directory)
