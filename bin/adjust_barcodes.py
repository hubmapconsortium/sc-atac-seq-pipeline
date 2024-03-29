#!/usr/bin/env pypy3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable

import add_barcodes_to_reads
import sciseq_add_barcodes_to_read_ids
from utils import Assay

ADJ_OUTPUT_DIR = Path("adj_fastq")
OUTPUT_FILENAME_PREFIX = "barcode_added"

adj_func_default = add_barcodes_to_reads.main
adj_funcs_special = {
    Assay.SCISEQ: sciseq_add_barcodes_to_read_ids.main,
}


def main(assay: Assay, input_dirs: Iterable[Path], output_filename_prefix, output_dir):
    ADJ_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    func = adj_funcs_special.get(assay, adj_func_default)
    print("Calling function", func, "to add barcodes")
    func(assay, input_dirs, output_filename_prefix, output_dir)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("directory", type=Path, nargs="+")
    p.add_argument("output_filename_prefix", nargs="?", default=OUTPUT_FILENAME_PREFIX)
    p.add_argument("output_dir", type=Path, nargs="?", default=ADJ_OUTPUT_DIR)
    args = p.parse_args()

    main(args.assay, args.directory, args.output_filename_prefix, args.output_dir)
