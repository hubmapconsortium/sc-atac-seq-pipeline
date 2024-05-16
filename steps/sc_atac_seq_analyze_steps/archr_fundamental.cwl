#!/usr/bin/env cwl-runner

class: CommandLineTool
id: archr_analyze
label: ArchR analysis
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38
  NetworkAccess:
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
      - $(inputs.bam_file)
      - $(inputs.bam_index)

inputs:
  bam_file:
    type: File
    inputBinding:
      position: 1
      prefix: --bam_file
      valueFrom: $(self.basename)
    doc: "The sorted BAM file with cell ids in the CB tag."

  bam_index:
    type: File
    inputBinding:
      position: 2
      prefix: --bam_index
      valueFrom: $(self.basename)
    doc: "The BAM file index."

  threads:
    type: int?
    inputBinding:
      position: 3
      prefix: --threads
    default: 2
    doc: "Number of threads to use."

outputs:
  cell_column_data_csv:
    type: File
    outputBinding:
      glob: "cell_column_data.csv"

  gene_row_data_csv:
    type: File
    outputBinding:
      glob: "gene_row_data.csv"

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
      glob: "cell_by_bin.mtx"

  cell_by_bin_barcodes:
    type: File
    outputBinding:
      glob: "barcodes.txt"

  cell_by_bin_bins:
    type: File
    outputBinding:
      glob: "bins.txt"

  rdata_file:
    type: File
    outputBinding:
      glob: "atacSeqStep1.RData"

baseCommand: [Rscript, /opt/run_ArchR_analysis_pt1.R]
