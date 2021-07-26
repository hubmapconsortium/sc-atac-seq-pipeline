#!/usr/bin/env python3
import json
import logging
from argparse import ArgumentParser
from collections import Counter
from os import fspath
from pathlib import Path
from statistics import median

import anndata
import HTSeq
import numpy as np

# In container
default_annotation_file = Path("/opt/supplementary-data/gencode.v32.annotation.gff3.gz")


def main(
    bam_file: Path,
    peak_file: Path,
    annotations_file: Path,
    cell_by_bin_file: Path = None,
):
    logging.info("Building exon/transcript index")
    # TODO: fix htseq                      ↓↓↓
    annotations_reader = HTSeq.GFF_Reader(fspath(annotations_file))
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
    transcripts = HTSeq.GenomicArrayOfSets("auto", stranded=True)

    for annotation in annotations_reader:
        if annotation.type == "exon":
            exons[annotation.iv] += annotation.attr["gene_id"]
        elif annotation.type == "transcript":
            transcripts[annotation.iv] += annotation.attr["gene_id"]

    logging.info("Building peak index")
    peaks = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    with open(peak_file) as f:
        for i, line in enumerate(f):
            name = f"peak_{i}"
            pieces = line.strip().split()
            chrom, start, stop = pieces[0], int(pieces[1]), int(pieces[2])
            iv = HTSeq.GenomicInterval(chrom, start, stop)
            peaks[iv] += name

    bam = HTSeq.BAM_Reader(fspath(bam_file))

    barcode_counts = Counter()
    barcode_reads_in_peaks = Counter()

    total_reads = 0
    mapped_reads = 0
    unmapped_reads = 0
    exonic_reads = 0
    intronic_reads = 0
    intergenic_reads = 0
    alignment_qualities = []

    for i, seg in enumerate(bam):
        if not (i % 10000):
            logging.debug(f"Processed {i} reads")
        barcode = seg.read.name.split(":", 1)[1]
        barcode_counts[barcode] += 1

        total_reads += 1
        if seg.aligned:
            mapped_reads += 1
            alignment_qualities.append(seg.aQual)
        else:
            unmapped_reads += 1
            continue

        read_exons = set()
        read_transcripts = set()

        for iv, val in exons[seg.iv].steps():
            read_exons |= val
        for iv, val in transcripts[seg.iv].steps():
            read_transcripts |= val

        read_introns = read_transcripts - read_exons

        if read_exons:
            exonic_reads += 1
        elif read_introns:
            intronic_reads += 1
        else:
            intergenic_reads += 1

        read_peaks = set()
        for iv, val in peaks[seg.iv].steps():
            read_peaks |= val
        barcode_reads_in_peaks[barcode] += bool(read_peaks)

    median_reads_in_peaks_mean = median(barcode_reads_in_peaks.values())

    proportion_mapped = mapped_reads / total_reads

    percentiles = [
        (0, "minimum"),
        (25, "25th_quantile"),
        (50, "median"),
        (75, "75th_quantile"),
        (100, "maximum"),
    ]
    five_number_dict = {key: np.percentile(alignment_qualities, pc) for pc, key in percentiles}

    qc_report = {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "unmapped_reads": unmapped_reads,
        "mapped_proportion": proportion_mapped,
        "mapping_quality": five_number_dict,
        "exonic_proportion": exonic_reads / total_reads,
        "intronic_proportion": intronic_reads / total_reads,
        "intergenic_proportion": intergenic_reads / total_reads,
        "median_proportion_reads_in_peaks": median_reads_in_peaks_mean,
    }

    if cell_by_bin_file:
        cell_by_bin = anndata.read_h5ad(cell_by_bin_file)
        qc_report["barcodes_passing_qc"] = cell_by_bin.shape[0]

    with open("qc_report.json", "w") as text_file:
        json.dump(qc_report, text_file, indent=4)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("bam_file", type=Path)
    p.add_argument("peak_file", type=Path)
    p.add_argument("cell_by_bin_file", type=Path, nargs="?")
    p.add_argument("--annotations-file", type=Path, default=default_annotation_file)
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format="%(asctime)-15s %(levelname)-7s | %(message)s")

    main(
        bam_file=args.bam_file,
        peak_file=args.peak_file,
        annotations_file=args.annotations_file,
        cell_by_bin_file=args.cell_by_bin_file,
    )
