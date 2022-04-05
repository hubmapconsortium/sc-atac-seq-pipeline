#!/usr/bin/env cwl-runner

class: CommandLineTool
id: archr_analyze
label: ArchR analysis
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.0b8
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

outputs:
  Fragment_Size_Distribution_pdf:
    type: File
    outputBinding:
      glob: "QualityControl/*/*-Fragment_Size_Distribution.pdf"

  TSS_by_Unique_Frags_pdf:
    type: File
    outputBinding:
      glob: "QualityControl/*/*-TSS_by_Unique_Frags.pdf"

  QC-Sample-FragSizes-TSSProfile_pdf:
    type: File
    outputBinding:
      glob: "*/Plots/QC-Sample-FragSizes-TSSProfile.pdf"

  QC-Sample-Statistics_pdf:
    type: File
    outputBinding:
      glob: "*/Plots/QC-Sample-Statistics.pdf"

  TSS-vs-Frags_pdf:
    type: File
    outputBinding:
      glob: "*/Plots/TSS-vs-Frags.pdf"

  Plot-UMAP-Sample-Clusters_pdf:
    type: File
    outputBinding:
      glob: "*/Plots/Plot-UMAP-Sample-Clusters.pdf"

  Peak-Call-Summary_pdf:
    type: File
    outputBinding:
      glob: "*/Plots/Peak-Call-Summary.pdf"

  GeneScores-Marker-Heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*/Plots/GeneScores-Marker-Heatmap.pdf"

  Peak-Marker-Heatmap_pdf:
    type: File
    outputBinding:
      glob: "*/Plots/Peak-Marker-Heatmap.pdf"

  peaks_csv:
    type: File
    outputBinding:
      glob: "peaks.csv"

  peak_markers_csv:
    type: File
    outputBinding:
      glob: "peak_markers.csv"

  peaks_bed:
    type: File
    outputBinding:
      glob: "peaks.bed"

  gene_markers_csv:
    type: File?
    outputBinding:
      glob: "gene_markers.csv"

  cell_column_data_csv:
    type: File
    outputBinding:
      glob: "cell_column_data.csv"

  gene_row_data_csv:
    type: File
    outputBinding:
      glob: "gene_row_data.csv"

  umap_coords_clusters_csv:
    type: File
    outputBinding:
      glob: "archr_umap_coords_clusters.csv"

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

baseCommand: [Rscript, /opt/run_ArchR_analysis.R]
