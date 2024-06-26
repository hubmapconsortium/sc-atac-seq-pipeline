#!/usr/bin/env cwl-runner

class: CommandLineTool
id: archr_init_analyze
label: ArchR initial analysis
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.1.4
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

  minTSS:
    type: double?
    inputBinding:
      position: 4
      prefix: --minTSS
    default: 1.5
    doc: "The minimum numeric transcription start site (TSS) enrichment score required to pass filtering. E.g. 1.5"

  minFrags:
    type: int?
    inputBinding:
      position: 5
      prefix: --minFrags
    default: 2000
    doc: "The minimum number of mapped ATAC-seq fragments required per cell to pass filtering. E.g. 2000"

  minCells:
    type: int?
    inputBinding:
      position: 6
      prefix: --minCells
    default: 1000
    doc: "The minimum number of cells in the ArchR project that must pass filtering before a warning message is printed. E.g. 1000"


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

  TSS_by_Unique_Frags_pdf:
    type: File
    outputBinding:
      glob: "QualityControl/*/*-TSS_by_Unique_Frags.pdf"

  Fragment_Size_Distribution_pdf:
    type: File
    outputBinding:
      glob: "QualityControl/*/*-Fragment_Size_Distribution.pdf"

  TSS-vs-Frags_pdf:
    type: File
    outputBinding:
      glob: "ArchRStep1/Plots/TSS-vs-Frags.pdf"

  QC-Sample-Statistics_pdf:
    type: File
    outputBinding:
      glob: "ArchRStep1/Plots/QC-Sample-Statistics.pdf"

  QC-Sample-FragSizes-TSSProfile_pdf:
    type: File
    outputBinding:
      glob: "ArchRStep1/Plots/QC-Sample-FragSizes-TSSProfile.pdf"

  image_file:
    type: File
    outputBinding:
      glob: "atacSeqStep1.RData"

  archr_project:
    type: Directory
    outputBinding:
      glob: "ArchRStep1"

baseCommand: [Rscript, /opt/run_ArchR_analysis_pt1.R]
