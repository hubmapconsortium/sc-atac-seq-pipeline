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
    umap_coords_csv: Path,
    cell_by_gene_raw_mtx: Path,
    cell_by_gene_smoothed_hdf5: Path,
    cell_by_bin_mtx: Path,
    cell_by_bin_barcodes: Path,
    cell_by_bin_bins: Path,
):
    umap_coords_df = pd.read_csv(umap_coords_csv, index_col=0)
    umap_coords_df.loc[:, "cluster"] = umap_coords_df.loc[:, "Clusters"].astype("category")
    umap_coords = umap_coords_df.loc[
        :, ["IterativeLSI#UMAP_Dimension_1", "IterativeLSI#UMAP_Dimension_2"]
    ].to_numpy()

    cell_by_bin_mat = scipy.io.mmread(cell_by_bin_mtx).astype(bool).tocsr()
    with open(cell_by_bin_barcodes) as f:
        barcodes = [line.strip() for line in f]
    with open(cell_by_bin_bins) as f:
        bins = [line.strip() for line in f]

    assert barcodes == list(umap_coords_df.index)

    obs = umap_coords_df.loc[:, ["cluster"]].copy()
    obsm = {"X_umap": umap_coords}

    chroms = []
    bin_start = []
    bin_stop = []
    for b in bins:
        chrom, pos = b.rsplit(":", 1)
        start, stop = (int(p) for p in pos.split("-"))
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

    cell_by_bin = anndata.AnnData(
        cell_by_bin_mat,
        obs=obs,
        obsm=obsm,
        var=var,
        dtype=bool,
    )

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
        obsm=obsm,
        var=pd.DataFrame(index=genes),
        layers={"smoothed": cell_by_gene_smoothed},
    )

    print("Saving cell by gene matrix")
    cell_by_gene.write_h5ad("cell_by_gene.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("umap_coords_csv", type=Path)
    p.add_argument("cell_by_gene_raw_mtx", type=Path)
    p.add_argument("cell_by_gene_smoothed_hdf5", type=Path)
    p.add_argument("cell_by_bin_mtx", type=Path)
    p.add_argument("cell_by_bin_barcodes", type=Path)
    p.add_argument("cell_by_bin_bins", type=Path)
    args = p.parse_args()

    main(
        umap_coords_csv=args.umap_coords_csv,
        cell_by_gene_raw_mtx=args.cell_by_gene_raw_mtx,
        cell_by_gene_smoothed_hdf5=args.cell_by_gene_smoothed_hdf5,
        cell_by_bin_mtx=args.cell_by_bin_mtx,
        cell_by_bin_barcodes=args.cell_by_bin_barcodes,
        cell_by_bin_bins=args.cell_by_bin_bins,
    )
