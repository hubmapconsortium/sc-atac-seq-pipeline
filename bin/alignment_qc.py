#!/usr/bin/env python3


import pysam
import json
from argparse import ArgumentParser
from pathlib import Path
import numpy as np

#Take a BAM file
#Read it with pysam
#Compute summary statistics of alignment quality
#And write those out to a JSON

def main(bam_file: Path, threads:int):

    bam = pysam.AlignmentFile(bam_file, "rb", threads=threads)

    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0
    alignment_qualities = []

    for seg in bam:
        total_reads += 1
        if seg.is_unmapped:
            unmapped_reads += 1
        else:
            mapped_reads += 1
            alignment_qualities.append(seg.mapping_quality)


    proportion_mapped = mapped_reads / total_reads

    percentiles = [0, 25, 50, 75, 100]
    sub_dict_keys = ["minimum", "25th_quantile", "median", "75th_quantile", "maximum"]
    five_number_dict = {sub_dict_keys[i]: np.percentile(alignment_qualities, percentiles[i]) for i in range(5)}

    qc_report = {"total_reads": total_reads, "mapped_reads": mapped_reads, "unmapped_reads": unmapped_reads, "proportion_mapped": proportion_mapped, \
    "mapping_quality": five_number_dict}

    with open("alignment_qc.json", "w") as text_file:
        json.dump(qc_report, text_file, indent=4)

if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('bam_file', type=Path)
    p.add_argument('threads', type=int)
    args = p.parse_args()

    main(args.bam_file, args.threads)
