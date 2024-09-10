#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse


def main(
    cell_column_data: Path,
    cell_by_gene_raw_mtx: Path,
    cell_by_gene_smoothed_hdf5: Path,
    cell_by_bin_mtx: Path,
    cell_by_bin_barcodes: Path,
    cell_by_bin_bins: Path,
    bin_size: int,
):

    cell_by_bin_mat = scipy.io.mmread(cell_by_bin_mtx).astype(bool).tocsr()
    # Barcode file looks like this:
    # BAM_data#GTCTGTCAAATCCGTCATGCCTAA
    # no header and assay name prefixed on barcode
    with open(cell_by_bin_barcodes) as f:
        barcodes = [line.strip() for line in f]
    # bins file looks like:
    # chr1 1 0
    # chr1 2 500
    with open(cell_by_bin_bins) as f:
        bins = [line.strip() for line in f]

    obs = pd.read_csv(cell_column_data, index_col=0)

    chroms = []
    bin_start = []
    bin_stop = []
    for b in bins:
        # chrom, pos = b.rsplit(":", 1)
        chrom, index, start = b.split()
        # Remove 'chr' from chromosome string
        prefix = "chr"
        if chrom.startswith(prefix):
            chrom = chrom[len(prefix) :]
        # print("chrom:{}, start:{}, class start:{}".format(chrom, start, type(start)))
        # start, stop = (int(p) for p in pos.split("-"))
        # 'start' is read to a string in expoential notation, so use float to convert it as int will fail
        # https://stackoverflow.com/questions/32861429/converting-number-in-scientific-notation-to-int
        stop = int(float(start)) + (
            bin_size - 1
        )  # E.g Archr is using bins of size 500 starting at 0
        chroms.append(chrom)
        bin_start.append(start)
        bin_stop.append(stop)
    var = pd.DataFrame(
        {
            "chrom": chroms,
            "bin_start": bin_start,
            "bin_stop": bin_stop,
        },
        index=bins,
    )

    print("obs:{}".format(obs))
    print("var:{}".format(var))

    cell_by_bin = anndata.AnnData(
        cell_by_bin_mat,
        obs=obs,
        var=var,
        dtype=bool,
    )
    cell_by_bin.obs = cell_by_bin.obs.dropna()
    print("Saving cell by bin matrix")
    cell_by_bin.write_h5ad("cell_by_bin.h5ad")

    cell_by_gene_raw = scipy.sparse.csr_matrix(scipy.io.mmread(cell_by_gene_raw_mtx))

    with h5py.File(cell_by_gene_smoothed_hdf5, "r") as f:
        cell_by_gene_smoothed = np.array(f["cell_by_gene_smoothed"]).T
        cells = [row.decode("utf-8") for row in np.array(f["barcodes"])]
        genes = [col.decode("utf-8") for col in np.array(f["genes"])]

    assert barcodes == cells

    cell_by_gene = anndata.AnnData(
        cell_by_gene_raw,
        obs=obs,
        var=pd.DataFrame(index=genes),
        layers={"smoothed": cell_by_gene_smoothed},
    )
    cell_by_gene.obs = cell_by_gene.obs.dropna()
    print("Saving cell by gene matrix")
    cell_by_gene.write_h5ad("cell_by_gene.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("cell_column_data", type=Path)
    p.add_argument("cell_by_gene_raw_mtx", type=Path)
    p.add_argument("cell_by_gene_smoothed_hdf5", type=Path)
    p.add_argument("cell_by_bin_mtx", type=Path)
    p.add_argument("cell_by_bin_barcodes", type=Path)
    p.add_argument("cell_by_bin_bins", type=Path)
    p.add_argument("bin_size", type=int)
    args = p.parse_args()

    main(
        cell_column_data=args.cell_column_data,
        cell_by_gene_raw_mtx=args.cell_by_gene_raw_mtx,
        cell_by_gene_smoothed_hdf5=args.cell_by_gene_smoothed_hdf5,
        cell_by_bin_mtx=args.cell_by_bin_mtx,
        cell_by_bin_barcodes=args.cell_by_bin_barcodes,
        cell_by_bin_bins=args.cell_by_bin_bins,
        bin_size=args.bin_size,
    )
