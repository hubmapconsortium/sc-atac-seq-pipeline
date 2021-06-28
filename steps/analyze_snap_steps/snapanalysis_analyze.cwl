#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_analyze
label: snap analysis analyze
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.1
  NetworkAccess:
    networkAccess: true

inputs:
  input_snap:
    type: File
    inputBinding:
      position: 1
      prefix: --input_snap
    doc: The selected barcodes object file.

  selected_barcodes:
    type: File?
    inputBinding:
      position: 2
      prefix: --selected_barcodes
    doc: The selected barcodes object file.

  encode_blacklist:
    type: File?
    inputBinding:
      position: 3
      prefix: --encode_blacklist
    doc: A BED file of ENCODE blacklist to prevent potential artifacts.

  gene_track:
    type: File?
    inputBinding:
      position: 4
      prefix: --gene_track
    doc: A BED file of gene tracks.

  gene_annotation:
    type: File?
    inputBinding:
      position: 5
      prefix: --gene_annotation
    doc: A GTF file of gene annotations.

  promoters:
    type: File?
    inputBinding:
      position: 6
      prefix: --promoters
    doc: A BED file of promoters.

  processes:
    type: int?
    inputBinding:
      position: 7
      prefix: --processes
    default: 1
    doc: Number of processes to use

outputs:
  snap_rds:
    type: File
    outputBinding:
      glob: "peaks_snap.rds"

  peaks_combined_bed:
    type: File
    outputBinding:
      glob: "peaks.combined.bed"

  CSV_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.csv"

  PDF_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.pdf"

  RDS_objects:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.rds"

  peaks_bed_file:
    type: File
    outputBinding:
      glob: "peaks.combined.bed"

  umap_coords_csv:
    type: File
    outputBinding:
      glob: "umap_coords_clusters.csv"

  cell_by_gene_raw_mtx:
    type: File
    outputBinding:
      glob: "cell_by_gene_raw.mtx"

  cell_by_gene_smoothed_hdf5:
    type: File
    outputBinding:
      glob: "cell_by_gene_smoothed.hdf5"

  cell_by_bin_mtx:
    type: File
    outputBinding:
      glob: "filtered_cell_by_bin.mtx"

  cell_by_bin_barcodes:
    type: File
    outputBinding:
      glob: "barcodes.txt"

  cell_by_bin_bins:
    type: File
    outputBinding:
      glob: "bins.txt"

baseCommand: [Rscript, /opt/snapAnalysis.R]
