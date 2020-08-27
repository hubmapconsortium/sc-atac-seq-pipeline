#!/usr/bin/env python3
from argparse import ArgumentParser
import json
from os import fspath
from pathlib import Path
from pprint import pprint
from typing import Iterable, Tuple

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
]

FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*_R1*.{extension}'
    for extension in FASTQ_EXTENSIONS:
        yield from directory.glob(pattern.format(extension=extension))

def find_bulk_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path]]:
    """
    :param directory:
    :return: Iterable of 2-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        r2_fastq_filename = r1_fastq_file.name.replace('_R1', '_R2')
        r2_fastq_file = r1_fastq_file.with_name(r2_fastq_filename)

        if r2_fastq_file.is_file():
            print(FOUND_PAIR_COLOR + 'Found two FASTQ files:' + NO_COLOR)
            print('\t', r1_fastq_file, sep='')
            print('\t', r2_fastq_file, sep='')
            yield r1_fastq_file, r2_fastq_file
        else:
            print(UNPAIRED_COLOR + 'Found unpaired FASTQ file:' + NO_COLOR)
            print('\t', r1_fastq_file, sep='')

def main(directories: Iterable[Path]):
    sequence_file_bundles = []

    for directory in directories:
        for r1_fastq_file, r2_fastq_file in find_bulk_fastq_files(directory):
            sequence_bundle = {
                'fastq_r1': {
                    'class': 'File',
                    'path': fspath(r1_fastq_file),
                },
                'fastq_r2': {
                    'class': 'File',
                    'path': fspath(r2_fastq_file),
                },
            }
            sequence_file_bundles.append(sequence_bundle)

    print('Sequence file bundles:')
    pprint(sequence_file_bundles)

    with open("input.json", "w") as text_file:
        json.dump(sequence_file_bundles, text_file)


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path, nargs='+')
    args = p.parse_args()

    main(args.directory)
