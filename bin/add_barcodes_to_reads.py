#!/usr/bin/env python3
import csv
import re
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path
from typing import Iterable, Mapping, Optional, Set

import barcodeutils as bu
from fastq_utils import Read, fastq_reader, find_grouped_fastq_files, revcomp

from utils import Assay

SNARESEQ_BARCODE_LENGTH = 8
SNARESEQ_BARCODE_STARTS = [0, 38, 76]
SNARESEQ_BARCODE_SEGMENTS = [
    slice(start, start + SNARESEQ_BARCODE_LENGTH) for start in SNARESEQ_BARCODE_STARTS
]

SNSEQ_BARCODE_SEGMENT = slice(0, 16)
# default for 10x, if not read from metadata
MULTIOME_10X_BARCODE_SEGMENT = slice(0, 16)

UMI_START = 84
UMI_LENGTH = 10
UMI_SEGMENT = slice(UMI_START, UMI_START + UMI_LENGTH)

metadata_filename_pattern = re.compile(r"^[0-9A-Fa-f]{32}-metadata.tsv$")


def create_and_output_new_read(orig_read: Read, new_read_id_prefix, output_file_obj):
    # Remove white space at beginning of ID (there shouldn't be any whitespace)
    orig_read_id_wo_left_whitespace = orig_read.read_id.lstrip()
    # Then remove the @ sign; it should be the first character and it will be put back later
    orig_read_id_wo_lw_or_ampersand = orig_read_id_wo_left_whitespace[1:]
    new_read_id = "@" + new_read_id_prefix + orig_read_id_wo_lw_or_ampersand
    new_read = Read(
        read_id=new_read_id,
        seq=orig_read.seq,
        unused=orig_read.unused,
        qual=orig_read.qual,
    )

    print(new_read.serialize(), file=output_file_obj)


def create_and_output_barcode_adjusted_reads(
    f1r: Read, f2r: Read, umi_seq, barcode_pieces, barf1addedout, barf2addedout
):
    new_read_id_prefix = "".join(barcode_pieces + [":", umi_seq, ":"])
    create_and_output_new_read(f1r, new_read_id_prefix, barf1addedout)
    create_and_output_new_read(f2r, new_read_id_prefix, barf2addedout)


def read_barcode_allowlist(barcode_filename: Path) -> Set[str]:
    with open(barcode_filename) as f:
        return set(f.read().split())


def find_metadata_file(directory: Path) -> Optional[Path]:
    """
    Finds and returns the first metadata file for a HuBMAP data set.
    Does not check whether the dataset ID (32 hex characters) matches
    the directory name, nor whether there might be multiple metadata files.
    """
    for file_path in directory.iterdir():
        if metadata_filename_pattern.match(file_path.name):
            return file_path


def main(
    assay: Assay,
    fastq_dirs: Iterable[Path],
    orig_fastq_dir: Path,
    output_filename_prefix,
    output_dir: Path,
):
    baraddedf1 = output_dir / f"{output_filename_prefix}_R1.fastq"
    baraddedf2 = output_dir / f"{output_filename_prefix}_R2.fastq"

    count = 0
    non_count = 0

    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, assay.fastq_count) for fastq_dir in fastq_dirs
    )

    metadata_file = find_metadata_file(orig_fastq_dir)
    if metadata_file is None:
        print("no metadata file found")
    else:
        print("metadata file found")

    with open(baraddedf1, "w") as barf1addedout, open(baraddedf2, "w") as barf2addedout:
        if assay == Assay.MULTIOME_10X:
            multiome_seg = MULTIOME_10X_BARCODE_SEGMENT
            barcode_filename = "/opt/atac_barcodes_rev.txt"
            barcode_allowlist = read_barcode_allowlist(barcode_filename)
            correcter = bu.BarcodeCorrecter(barcode_allowlist, edit_distance=1)
            if metadata_file is not None and metadata_file.is_file():
                with open(metadata_file, newline="") as f:
                    r = csv.DictReader(f, delimiter="\t")
                    metadata = next(r)
                    # defaults
                    offset = 0
                    length = 16
                    if (
                        "cell_barcode_offset" in metadata
                        and metadata["cell_barcode_offset"].isdigit()
                    ):
                        offset = int(metadata["cell_barcode_offset"])
                        print(f"offset is {offset}")
                    if "cell_barcode_size" in metadata and metadata["cell_barcode_size"].isdigit():
                        length = int(metadata["cell_barcode_size"])
                        print(f"length is {length}")
                    multiome_seg = slice(offset, offset + length)

        for fastq1_file, fastq2_file, barcode_file in all_fastqs:
            i = 0
            print("Adding barcodes to", fastq1_file, "and", fastq2_file, "using", barcode_file)
            fastq1_reader = fastq_reader(fastq1_file)
            fastq2_reader = fastq_reader(fastq2_file)
            barcode_reader = fastq_reader(barcode_file)

            for i, (f1r, f2r, bar) in enumerate(
                zip(fastq1_reader, fastq2_reader, barcode_reader), 1
            ):
                if assay == Assay.SNARESEQ:
                    # This code assumes the barcode fastq is the third fastq; e.g. '...R3...'
                    barcode_pieces = [revcomp(bar.seq[s]) for s in SNARESEQ_BARCODE_SEGMENTS]
                    # For some reason the barcode pieces are reversed in the list
                    # so get them back in the right order
                    barcode_pieces.reverse()
                    umi_seq = bar.seq[UMI_SEGMENT]
                elif assay == Assay.SNSEQ:
                    barcode_pieces = [bar.seq[SNSEQ_BARCODE_SEGMENT]]
                    umi_seq = ""
                elif assay == Assay.MULTIOME_10X:
                    barcode_pieces = [bar.seq[multiome_seg]]
                    barcode_pieces = [correcter.correct(barcode) for barcode in barcode_pieces]

                    if barcode_pieces[0] is None:
                        barcode_pieces = ["N" * 16]
                        non_count += 1
                    else:
                        count += 1

                    umi_seq = ""
                else:
                    print("Could not adjust barcodes for assay {}".format(assay))

                create_and_output_barcode_adjusted_reads(
                    f1r,
                    f2r,
                    umi_seq,
                    barcode_pieces,
                    barf1addedout,
                    barf2addedout,
                )
    if assay == Assay.MULTIOME_10X:
        print(f"added {count} barcodes, {non_count} barcodes did not match")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", type=Assay)
    p.add_argument("fastq_dirs", type=Path, nargs="+")
    p.add_argument("orig_fastq_dir", type=Path)
    p.add_argument("output_filename_prefix")
    p.add_argument("output_dir", type=Path)

    args = p.parse_args()

    main(
        args.assay,
        args.fastq_dirs,
        args.orig_fastq_dir,
        args.output_filename_prefix,
        args.output_dir,
    )
