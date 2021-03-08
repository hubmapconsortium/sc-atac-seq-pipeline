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
    cell_by_gene_matrix: Path,
    cell_by_bin_mtx: Path,
    cell_by_bin_barcodes: Path,
    cell_by_bin_bins: Path,
):
    umap_coords_df = pd.read_csv(umap_coords_csv, index_col=0)
    umap_coords_df.loc[:, "cluster"] = umap_coords_df.loc[:, "cluster"].astype("category")
    umap_coords = umap_coords_df.loc[:, ["umap.1", "umap.2"]].to_numpy()

    cell_by_bin_mat = scipy.io.mmread(cell_by_bin_mtx).astype(bool).tocsr()
    with open(cell_by_bin_barcodes) as f:
        barcodes = [line.strip() for line in f]
    with open(cell_by_bin_bins) as f:
        bins = [line.strip() for line in f]

    assert barcodes == list(umap_coords_df.index)

    obs = umap_coords_df.loc[:, ["cluster"]].copy()
    obsm = {"X_umap": umap_coords}

    cell_by_bin = anndata.AnnData(
        cell_by_bin_mat,
        obs=obs,
        obsm=obsm,
        var=pd.DataFrame(index=bins),
        dtype=bool,
    )

    print("Saving cell by bin matrix")
    cell_by_bin.write_h5ad("cell_by_bin.h5ad")

    with h5py.File(cell_by_gene_matrix, "r") as f:
        cell_by_gene_raw = scipy.sparse.csr_matrix(f["cell_by_gene_raw"])
        cell_by_gene_smoothed = np.array(f["cell_by_gene"])
        cells = [row.decode("utf-8") for row in np.array(f["row_names"])]
        genes = [col.decode("utf-8") for col in np.array(f["col_names"])]

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
    p.add_argument("cell_by_gene_matrix", type=Path)
    p.add_argument("cell_by_bin_mtx", type=Path)
    p.add_argument("cell_by_bin_barcodes", type=Path)
    p.add_argument("cell_by_bin_bins", type=Path)
    args = p.parse_args()
