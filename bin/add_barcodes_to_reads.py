#!/usr/bin/env python3
from argparse import ArgumentParser
from itertools import chain
from pathlib import Path
from typing import Iterable, Mapping, Set

from fastq_utils import Read, fastq_reader, find_grouped_fastq_files, revcomp

BARCODE_LENGTH = 8
BARCODE_STARTS = [0, 38, 76]
BARCODE_SEGMENTS = [slice(start, start + BARCODE_LENGTH) for start in BARCODE_STARTS]
UMI_START = 84
UMI_LENGTH = 10
UMI_SEGMENT = slice(UMI_START, UMI_START + UMI_LENGTH)

def main(
    fastq_dirs: Iterable[Path],
    output_dir: Path = Path()
):

    baraddedf1 = output_dir / "barcode_added_R1.fastq"
    baraddedf2 = output_dir / "barcode_added_R2.fastq"

    all_fastqs = chain.from_iterable(
        find_grouped_fastq_files(fastq_dir, 3) for fastq_dir in fastq_dirs
    )

    with open(baraddedf1, "w") as barf1addedout, open(baraddedf2, "w") as barf2addedout:
        # This code assumes the barcode fastq is the third fastq; e.g. '...R3...'
        for fastq1_file, fastq2_file, barcode_file in all_fastqs:
            usable_count = 0
            i = 0
            print("Adding barcodes to", fastq1_file, "and", fastq2_file, "using", barcode_file)
            fastq1_reader = fastq_reader(fastq1_file)
            fastq2_reader = fastq_reader(fastq2_file)
            barcode_reader = fastq_reader(barcode_file)

            for i, (f1r, f2r, bar) in enumerate(zip(fastq1_reader, fastq2_reader, barcode_reader), 1):
                barcode_pieces = [revcomp(bar.seq[s]) for s in BARCODE_SEGMENTS]
                # For some reason the barcoe pieces are reversed in the list
                # so get them back in the right order
                barcode_pieces.reverse()
                #print("barcode pieces={}".format(barcode_pieces))
                
                umi_seq = bar.seq[UMI_SEGMENT]
                new_read_id_prefix = "".join(barcode_pieces + [":", umi_seq, ":"])

                # Remove white space at beginning of ID (there shouldn't be any whitespace)
                old_f1r_id_wo_left_whitespace = f1r.read_id.lstrip()
                # Then remove the @ sign; it should be the first character and it will be put back later
                old_f1r_id_wo_lw_or_ampersand = old_f1r_id_wo_left_whitespace[1:]
                new_f1r_id = '@' + new_read_id_prefix + old_f1r_id_wo_lw_or_ampersand
                new_f1r = Read(
                        read_id=new_f1r_id,
                        seq=f1r.seq,
                        unused=f1r.unused,
                        qual=f1r.qual,
                    )

                # Remove white space at beginning of ID (there shouldn't be any whitespace)
                old_f2r_id_wo_left_whitespace = f2r.read_id.lstrip()
                # Then remove the @ sign; it should be the first character and it will be put back later
                old_f2r_id_wo_lw_or_ampersand = old_f2r_id_wo_left_whitespace[1:]
                new_f2r_id = '@' + new_read_id_prefix + old_f2r_id_wo_lw_or_ampersand
                new_f2r = Read(
                        read_id=new_f2r_id,
                        seq=f2r.seq,
                        unused=f2r.unused,
                        qual=f2r.qual,
                    )

                print(new_f1r.serialize(), file=barf1addedout)
                print(new_f2r.serialize(), file=barf2addedout)

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("fastq_dirs", type=Path, nargs="+")
    args = p.parse_args()

    main(
        fastq_dirs=args.fastq_dirs
    )

