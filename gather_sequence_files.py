#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from os import fspath
from typing import Iterable, Tuple
from string import Template
import json

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
]



SEQUENCES_TEMPLATE = Template("""{
    "input_fastq1": {
      "class": "File",
      "path": "$fastq1"
    },
    "input_fastq2": {
      "class": "File",
      "path": "$fastq2"
    },
    "input_barcode_fastq": {
      "class": "File",
      "path": "$fastq3"
    }
}""")

BULK_SEQUENCES_TEMPLATE = Template("""{
    "input_fastq1": {
      "class": "File",
      "path": "$fastq1"
    },
    "input_fastq2": {
      "class": "File",
      "path": "$fastq2"
    },
    "input_barcode_fastq": {
      "class": "File",
      "path": "$fastq3"
    }
}""")


FOUND_PAIR_COLOR = '\033[01;32m'
UNPAIRED_COLOR = '\033[01;31m'
NO_COLOR = '\033[00m'

def find_r1_fastq_files(directory: Path) -> Iterable[Path]:
    pattern = '**/*_R1*.{extension}'
    print("in find_r1_fastq_files with path:{}".format(directory))
    for extension in FASTQ_EXTENSIONS:
        print("yielding files from directory {}".format(directory))
        yield from directory.glob(pattern.format(extension=extension))

def find_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path, Path]]:
    """
    Specific to 10X FASTQ filename conventions. Returns all paired R1/R2
    FASTQ files in any subdirectory of 'directory'.

    :param directory:
    :return: Iterable of 3-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
     [2] R3 FASTQ file
"""
    for r1_fastq_file in find_r1_fastq_files(directory):
        r2_fastq_filename = r1_fastq_file.name.replace('_R1', '_R2')
        r3_fastq_filename = r1_fastq_file.name.replace('_R1', '_R3')
        r2_fastq_file = r1_fastq_file.with_name(r2_fastq_filename)
        r3_fastq_file = r1_fastq_file.with_name(r3_fastq_filename)

        if r2_fastq_file.is_file() and r3_fastq_file.is_file():
            print(FOUND_PAIR_COLOR + 'Found three FASTQ files:' + NO_COLOR)
            print('\t', r1_fastq_file, sep='')
            print('\t', r2_fastq_file, sep='')
            print('\t', r3_fastq_file, sep='')
            yield r1_fastq_file, r2_fastq_file, r3_fastq_file
        else:
            print(UNPAIRED_COLOR + 'Found unpaired FASTQ file:' + NO_COLOR)
            print('\t', r1_fastq_file, sep='')


def find_bulk_fastq_files(directory: Path) -> Iterable[Tuple[Path, Path]]:
    """
    Specific to 10X FASTQ filename conventions. Returns all paired R1/R2
    FASTQ files in any subdirectory of 'directory'.

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



def main(directory: Path, bulk: bool):
    sequence_file_bundles = []

    if bulk:
        for r1_fastq_file, r2_fastq_file in find_bulk_fastq_files(directory):
            r1_fastq_path = fspath(r1_fastq_file)
            r2_fastq_path = fspath(r2_fastq_file)

            json_template_dict = { 'fastq1':r1_fastq_path, 'fastq2':r2_fastq_path}
            json_template_with_substitutes = BULK_SEQUENCES_TEMPLATE.substitute(json_template_dict)

            # Convert string representation of dictionary to an python dictionary
            # https://stackoverflow.com/questions/988228/convert-a-string-representation-of-a-dictionary-to-a-dictionary
            sequence_bundle = json.loads(json_template_with_substitutes)

            #print("sequence bundle is:", json.dumps(sequence_bundle))
            sequence_file_bundles.append(sequence_bundle)
        with open("input.json", "w") as text_file:
            json.dump(sequence_file_bundles, text_file)

    else:
        for r1_fastq_file, r2_fastq_file, r3_fastq_file in find_fastq_files(directory):

            r1_fastq_path = fspath(r1_fastq_file)
            r2_fastq_path = fspath(r2_fastq_file)
            r3_fastq_path = fspath(r3_fastq_file)

            json_template_dict = { 'fastq1':r1_fastq_path, 'fastq2':r2_fastq_path, 'fastq3':r3_fastq_path }
            json_template_with_substitutes = SEQUENCES_TEMPLATE.substitute(json_template_dict)

            # Convert string representation of dictionary to an python dictionary
            # https://stackoverflow.com/questions/988228/convert-a-string-representation-of-a-dictionary-to-a-dictionary
            sequence_bundle = json.loads(json_template_with_substitutes)

            #print("sequence bundle is:", json.dumps(sequence_bundle))
            sequence_file_bundles.append(sequence_bundle)
        with open("input.json", "w") as text_file:
            json.dump(sequence_file_bundles, text_file)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('directory', type=Path)
    p.add_argument('bulk', nargs='?', type=bool)
    args = p.parse_args()

    main(args.directory, main.bulk)
