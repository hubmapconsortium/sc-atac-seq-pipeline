#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path
from typing import Iterable, Mapping, Set

from utils import Assay
from fastq_utils import Read, fastq_reader, find_grouped_fastq_files, revcomp

SNARESEQ_BARCODE_LENGTH = 8
SNARESEQ_BARCODE_STARTS = [0, 38, 76]
SNARESEQ_BARCODE_SEGMENTS = [slice(start, start + SNARESEQ_BARCODE_LENGTH) for start in SNARESEQ_BARCODE_STARTS]

SCISEQ_BARCODE_SEGMENT = slice(0, 16)

UMI_START = 84
UMI_LENGTH = 10
UMI_SEGMENT = slice(UMI_START, UMI_START + UMI_LENGTH)

def create_and_output_new_read(orig_read: Read, new_read_id_prefix, output_file_obj):
    # Remove white space at beginning of ID (there shouldn't be any whitespace)
    orig_read_id_wo_left_whitespace = orig_read.read_id.lstrip()
    # Then remove the @ sign; it should be the first character and it will be put back later
    orig_read_id_wo_lw_or_ampersand = orig_read_id_wo_left_whitespace[1:]
    new_read_id = '@' + new_read_id_prefix + orig_read_id_wo_lw_or_ampersand
    new_read = Read(
            read_id=new_read_id,
            seq=orig_read.seq,
            unused=orig_read.unused,
            qual=orig_read.qual,
        )

    print(new_read.serialize(), file=output_file_obj)

def create_and_output_barcode_adjusted_reads(f1r: Read, f2r: Read, umi_seq, barcode_pieces, barf1addedout, barf2addedout):
    new_read_id_prefix = "".join(barcode_pieces + [":", umi_seq, ":"])
    create_and_output_new_read(f1r, new_read_id_prefix, barf1addedout)
    create_and_output_new_read(f2r, new_read_id_prefix, barf2addedout)


def main(
    assay: Assay,
    fastq_dirs: Iterable[Path],
    output_dir: Path = Path()
):

    baraddedf1 = output_dir / "barcode_added_R1.fastq"
    baraddedf2 = output_dir / "barcode_added_R2.fastq"

    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, assay.fastq_count) for fastq_dir in fastq_dirs
    )

    with open(baraddedf1, "w") as barf1addedout, open(baraddedf2, "w") as barf2addedout:
            #for fastq1_file, fastq2_file, barcode_file in all_fastqs:
            for fastq1_file, barcode_file, fastq2_file, in all_fastqs:

                i = 0
                print("Adding barcodes to", fastq1_file, "and", fastq2_file, "using", barcode_file)
                fastq1_reader = fastq_reader(fastq1_file)
                fastq2_reader = fastq_reader(fastq2_file)
                barcode_reader = fastq_reader(barcode_file)
    
                for i, (f1r, f2r, bar) in enumerate(zip(fastq1_reader, fastq2_reader, barcode_reader), 1):
                    if (assay == Assay.SNARESEQ):
                        # This code assumes the barcode fastq is the third fastq; e.g. '...R3...'
                        barcode_pieces = [revcomp(bar.seq[s]) for s in SNARESEQ_BARCODE_SEGMENTS]
                        # For some reason the barcode pieces are reversed in the list
                        # so get them back in the right order
                        barcode_pieces.reverse()
                        umi_seq = bar.seq[UMI_SEGMENT]
                    elif (assay == Assay.SNSEQ):
                        barcode_pieces = [f2r.seq[SCISEQ_BARCODE_SEGMENT]]
                        #print("barcode pieces={}".format(barcode_pieces))
                        umi_seq = ""
                    else:
                        print("Could not adjust barcodes for assay {}".format(assay))

                    create_and_output_barcode_adjusted_reads(f1r, f2r, umi_seq, barcode_pieces, barf1addedout, barf2addedout)

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("assay", type=Assay)
    p.add_argument("fastq_dirs", type=Path, nargs="+")

    args = p.parse_args()

    main(
        fastq_dirs=args.fastq_dirs
    )

