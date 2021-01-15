#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from subprocess import run
from typing import Iterable

from fastq_utils import find_grouped_fastq_files

from utils import Assay

CONCAT_OUTPUT_DIR = Path("concat_output_dir")

# No point in involving the Python interpreter in this decompression;
# it'll be a lot faster just to pipe 'gunzip -c' to each file
GUNZIP_COMMAND_TEMPLATE = [
    "gunzip",
    "-c",
    "{fastq_file}",
]
MERGED_FASTQ_R1 = Path(str(CONCAT_OUTPUT_DIR) + "/merged_R1.fastq")
MERGED_FASTQ_R2 = Path(str(CONCAT_OUTPUT_DIR) + "/merged_R2.fastq")
MERGED_FASTQ_BARCODE = Path(str(CONCAT_OUTPUT_DIR) + "/merged_R3.fastq")


def decompress_concat_fastq(input_fastq: Path, merged_fastq: Path):
    print("Concatenating", input_fastq, "to", merged_fastq)
    with open(merged_fastq, "ab") as out:
        command = [piece.format(fastq_file=input_fastq) for piece in GUNZIP_COMMAND_TEMPLATE]
        run(command, stdout=out, check=True)


def main(directories: Iterable[Path], assay: Assay):
    CONCAT_OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
    for directory in directories:
        for fastqs_list in find_grouped_fastq_files(directory, assay.fastq_count):
            num_fastqs = len(fastqs_list)
            if num_fastqs == 3:
                # For assays that use three files
                # flipping R2 and R3 here is deliberate. For the input data, barcodes
                # are in R2, and R1 and R3 are the two reads which should be aligned
                # to the genome. Internally in this pipeline, it's more straightforward
                # to rename these: R1 -> R1, R2 -> barcode, R3 -> R2
                r1_fastq = fastqs_list[0]
                decompress_concat_fastq(r1_fastq, MERGED_FASTQ_R1)
                r2_fastq = fastqs_list[1]
                decompress_concat_fastq(r2_fastq, MERGED_FASTQ_BARCODE)
                r3_fastq = fastqs_list[2]
                decompress_concat_fastq(r3_fastq, MERGED_FASTQ_R2)
            elif num_fastqs == 2:
                r1_fastq = fastqs_list[0]
                decompress_concat_fastq(r1_fastq, MERGED_FASTQ_R1)
                r2_fastq = fastqs_list[1]
                decompress_concat_fastq(r2_fastq, MERGED_FASTQ_R2)
            else:
                print(
                    "Could not unzip and concatenate fastqs becuase there are {} of them".format(
                        num_fastqs
                    )
                )
                exit(1)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("directory", type=Path, nargs="+")
    p.add_argument("assay", type=Assay)
    args = p.parse_args()

    main(args.directory, args.assay)
