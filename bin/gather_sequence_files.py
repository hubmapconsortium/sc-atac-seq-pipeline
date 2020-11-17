#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from os import fspath
from pathlib import Path

from fastq_utils import find_grouped_fastq_files

from utils import Assay


def main(assay: Assay, directory: Path):
    sequence_file_bundles = []

    for fastq_files in find_grouped_fastq_files(directory, assay.fastq_count, verbose=True):
        sequence_bundle = {
            f"input_fastq{n}": {"class": "File", "path": fspath(fastq_file)}
            for n, fastq_file in enumerate(fastq_files, 1)
        }
        sequence_file_bundles.append(sequence_bundle)

    with open("input.json", "w") as text_file:
        json.dump(sequence_file_bundles, text_file)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("directory", type=Path)
    args = p.parse_args()

    main(args.assay, args.directory)
