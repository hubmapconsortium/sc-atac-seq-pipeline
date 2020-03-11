from argparse import ArgumentParser
import gzip
from pathlib import Path

UMI_LENGTH = 8
UMI_START = 16

BARCODE_LENGTH = 10
LIGATION_BARCODE_START = 0
RT_BARCODE_START = 24

# Specific format for SnapTools
READ_ID_FORMAT = '@{barcode}:{umi}:{previous_read_id}'

def attach_barcode_umi_read2_id(
        fastq_r1: Path,
        fastq_r2: Path,
        output_file: Path
):
    f1 = gzip.open(fastq_r1)
    f2 = gzip.open(fastq_r2)
    f3 = gzip.open(output_file, 'wb')

    line1 = f1.readline()
    line2 = f2.readline()
    total_line = 0
    filtered_line = 0

    while line1:
        line1 = f1.readline().lstrip('@')
        # print("read1: ", line1)
        # first check if the ligation barcode match with the barcode

        ligation_barcode = line1[LIGATION_BARCODE_START:LIGATION_BARCODE_START + BARCODE_LENGTH]
        rt_barcode = line1[RT_BARCODE_START:RT_BARCODE_START + BARCODE_LENGTH]
        umi = line1[UMI_START:UMI_START + UMI_LENGTH]

        new_read_id = READ_ID_FORMAT.format(
            barcode=ligation_barcode + rt_barcode,
            umi=umi,
            previous_read_id=line1,
        )
        f3.write(new_read_id)

        second_line = f2.readline()
        f3.write(second_line)

        third_line = f2.readline()
        f3.write(third_line)

        four_line = f2.readline()
        f3.write(four_line)

        line2 = f2.readline()

        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()

    f1.close()
    f2.close()
    f3.close()

if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument('fastq_r1', type=Path)
    p.add_argument('fastq_r2', type=Path)
    p.add_argument('output_file', type=Path)
    args = p.parse_args()

    attach_barcode_umi_read2_id(args.fastq_r1, args.fastq_r2, args.output_path)
