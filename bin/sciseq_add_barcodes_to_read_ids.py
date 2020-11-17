#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Union

from fastq_utils import Read, fastq_reader, smart_open

# Specific format for SnapTools
READ_ID_FORMAT = '@{barcode}:{umi}:{previous_read_id}'

def convert_fastq(input_path: Path, output_path: Union[str, Path]):
    print('Converting', input_path, 'to', output_path)
    with smart_open(output_path, 'wt') as f:
        for read in fastq_reader(input_path):
            # We know the read ID format of this single data set, but
            # be lenient in what we accept and how we parse this
            barcode, *rest = read.read_id.lstrip('@').split(':')
            new_read_id = READ_ID_FORMAT.format(
                barcode=barcode,
                umi='',
                previous_read_id=''.join(rest),
            )

            new_read = Read(
                read_id=new_read_id,
                seq=read.seq,
                unused=read.unused,
                qual=read.qual,
            )
            print(new_read.serialize(), file=f)

def main(fastq_r1: Path, fastq_r2: Path, output_filename_prefix: str):
    convert_fastq(fastq_r1, f'{output_filename_prefix}.R1.fastq')
    convert_fastq(fastq_r2, f'{output_filename_prefix}.R2.fastq')

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument('fastq_r1', type=Path)
    p.add_argument('fastq_r2', type=Path)
    p.add_argument('output_filename_prefix')
    args = p.parse_args()

    main(args.fastq_r1, args.fastq_r2, args.output_filename_prefix)
