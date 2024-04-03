#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path
from typing import Iterable, Union

from fastq_utils import Read, fastq_reader, find_grouped_fastq_files

from utils import Assay

READ_ID_FORMAT = "@{barcode}:{umi}:{previous_read_id}"


def convert_fastq(input_path: Path, output_path: Union[str, Path]):
    print("Converting", input_path, "to", output_path)
    with open(output_path, "wt") as f:
        for read in fastq_reader(input_path):
            # We know the read ID format of this single data set, but
            # be lenient in what we accept and how we parse this
            barcode, *rest = read.read_id.lstrip("@").split(":")
            new_read_id = READ_ID_FORMAT.format(
                barcode=barcode,
                umi="",
                previous_read_id="".join(rest),
            )

            new_read = Read(
                read_id=new_read_id,
                seq=read.seq,
                unused=read.unused,
                qual=read.qual,
            )
            print(new_read.serialize(), file=f)


def main(
    assay: Assay,
    fastq_dirs: Iterable[Path],
    orig_dir: Path,
    output_filename_prefix: str,
    output_dir: Path,
):
    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, assay.fastq_count) for fastq_dir in fastq_dirs
    )

    for fastq_r1, fastq_r2 in all_fastqs:
        convert_fastq(fastq_r1, output_dir / f"{output_filename_prefix}_R1.fastq")
        convert_fastq(fastq_r2, output_dir / f"{output_filename_prefix}_R2.fastq")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", type=Assay)
    p.add_argument("output_filename_prefix")
    p.add_argument("output_dir", type=Path)
    p.add_argument("fastq_dirs", type=Path, nargs="+")

    args = p.parse_args()
    main(args.assay, args.fastq_dirs, args.output_filename_prefix, args.output_dir)
