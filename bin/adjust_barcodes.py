#!/usr/bin/env pypy3
from argparse import ArgumentParser
from pathlib import Path
from typing import Iterable
from utils import Assay

ADJ_OUTPUT_DIR = Path("adj_fastq")

adj_funcs = {
    #Assay.SLIDESEQ: extract_slideseq_barcodes.main,
    Assay.SNSEQ: adjust_barcodes.main
}

def main(assay: Assay, input_dirs: Iterable[Path]):
    ADJ_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    if assay in adj_funcs:
        print("Calling function to add barcodes")
        adj_funcs[assay](input_dirs, output_dir=ADJ_OUTPUT_DIR)
    else:
        print("No barcode adjustment to perform for assay", assay)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", choices=list(Assay), type=Assay)
    p.add_argument("directory", type=Path, nargs="+")
    args = p.parse_args()

    main(args.assay, args.directory)

