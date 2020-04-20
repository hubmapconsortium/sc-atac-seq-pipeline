#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from os import fspath
from subprocess import check_call
from typing import Iterable, Tuple
from string import Template

FASTQ_EXTENSIONS = [
    'fastq',
    'fastq.gz',
]

# cwltool --debug --leave-tmpdir --tmpdir-prefix=/mnt/tmp/ --tmp-outdir-prefix=/mnt/tmp/ ../create_snap_and_analyze.cwl ../create_snap_and_analyze.openstack.human.json
SC_ATAC_SEQ_COMMAND = [
    'cwltool',
    '--tmpdir-prefix=/mnt/tmp/',
    '--tmp-outdir-prefix=/mnt/tmp/',
    '/mnt/gitroot/sc-atac-seq-pipeline/create_snap_and_analyze.cwl',
    'input.json'
]


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
    :return: Iterable of 2-tuples:
     [0] R1 FASTQ file
     [1] R2 FASTQ file
    """
    for r1_fastq_file in find_r1_fastq_files(directory):
        print("file is:{}".format(r1_fastq_file))
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

def main(threads: int, directory: Path):
    command = [
        piece.format(threads=threads)
        for piece in SC_ATAC_SEQ_COMMAND
    ]
    print("command is:{}".format(command))
    for r1_fastq_file, r2_fastq_file, r3_fastq_file in find_fastq_files(directory):

        with open("/mnt/gitroot/sc-atac-seq-pipeline/input_template_json.txt", "r") as input_json_template_file:
            input_json_str = input_json_template_file.read()
        print("json template is {}".format(input_json_str))
        json_template = Template(input_json_str)

        r1_fastq_path = fspath(r1_fastq_file)
        r2_fastq_path = fspath(r2_fastq_file)
        r3_fastq_path = fspath(r3_fastq_file)
        json_template_dict = { 'fastq1':r1_fastq_path, 'fastq2':r3_fastq_path, 'fastq3':r2_fastq_path }
        json_template_with_substitutes = json_template.safe_substitute(json_template_dict)
        with open("input.json", "w") as text_file:
                text_file.write("{}".format(json_template_with_substitutes))
        print("input json is:{}".format(json_template_with_substitutes))

        print("Running:{}".format(' '.join(command)))
        #check_call(command)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('-p', '--threads', type=int)
    p.add_argument('directory', type=Path)
    args = p.parse_args()

    main(args.threads, args.directory)
