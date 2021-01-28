#!/usr/bin/env pypy3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable

import add_barcodes_to_reads
import sciseq_add_barcodes_to_read_ids
from utils import Assay

ADJ_OUTPUT_DIR = Path("adj_fastq")
OUTPUT_FILENAME_PREFIX = "barcode_added"

adj_funcs = {
    Assay.SNARESEQ: add_barcodes_to_reads.main,
    Assay.SCISEQ: sciseq_add_barcodes_to_read_ids.main,
    Assay.SNSEQ: add_barcodes_to_reads.main,
}


def main(assay: Assay, input_dirs: Iterable[Path], output_filename_prefix, output_dir):
    ADJ_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    if assay in adj_funcs:
        print("Calling function to add barcodes")
        adj_funcs[assay](assay, input_dirs, output_filename_prefix, output_dir)
    else:
        print("No barcode adjustment to perform for assay", assay)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("directory", type=Path, nargs="+")
    p.add_argument("output_filename_prefix", nargs="?", default=OUTPUT_FILENAME_PREFIX)
    p.add_argument("output_dir", type=Path, nargs="?", default=ADJ_OUTPUT_DIR)
    args = p.parse_args()

    main(args.assay, args.directory, args.output_filename_prefix, args.output_dir)
