cwlVersion: v1.0
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:latest
baseCommand: /opt/convert_to_h5ad.py

inputs:
  cell_column_data:
    type: File
    inputBinding:
      position: 0
  cell_by_gene_raw_mtx:
    type: File
    inputBinding:
      position: 1
  cell_by_gene_smoothed_hdf5:
    type: File
    inputBinding:
      position: 2
  cell_by_bin_mtx:
    type: File
    inputBinding:
      position: 3
  cell_by_bin_barcodes:
    type: File
    inputBinding:
      position: 4
  cell_by_bin_bins:
    type: File
    inputBinding:
      position: 5

  bin_size:
    type: int?
    inputBinding:
      position: 6
    default: 500

outputs:
  cell_by_bin_h5ad:
    type: File
    outputBinding:
      glob: 'cell_by_bin.h5ad'
  cell_by_gene_h5ad:
    type: File
    outputBinding:
      glob: 'cell_by_gene.h5ad'

