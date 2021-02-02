#!/usr/bin/env python3
import json
from argparse import ArgumentParser
from os import fspath
from pathlib import Path
from pprint import pprint
from typing import Iterable

from fastq_utils import find_grouped_fastq_files


def main(directories: Iterable[Path]):
    sequence_file_bundles = []

    for directory in directories:
        for r1_fastq_file, r2_fastq_file in find_grouped_fastq_files(directory, 2):
            sequence_bundle = {
                "fastq_r1": {
                    "class": "File",
                    "path": fspath(r1_fastq_file),
                },
                "fastq_r2": {
                    "class": "File",
                    "path": fspath(r2_fastq_file),
                },
            }
            sequence_file_bundles.append(sequence_bundle)

    print("Sequence file bundles:")
    pprint(sequence_file_bundles)

    with open("input.json", "w") as text_file:
        json.dump(sequence_file_bundles, text_file)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path, nargs="+")
    args = p.parse_args()

    main(args.directory)
