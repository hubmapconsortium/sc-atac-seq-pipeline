#!/usr/bin/env cwl-runner
class: CommandLineTool
id: archr_clustering
label: ArchR clustering
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

  image_file:
    type: File
    inputBinding:
      position: 3
      prefix: --image
      valueFrom: $(self.basename)
    doc: "The R image from the previous archR analysis step"

  threads:
    type: int?
    inputBinding:
      position: 4
      prefix: --threads
    default: 2
    doc: "Number of threads to use."

  minTSS:
    type: double?
    inputBinding:
      position: 5
      prefix: --minTSS
    default: 1.5
    doc: "The minimum numeric transcription start site (TSS) enrichment score required to pass filtering. E.g. 1.5"

  minFrags:
    type: int?
    inputBinding:
      position: 6
      prefix: --minFrags
    default: 2000
    doc: "The minimum number of mapped ATAC-seq fragments required per cell to pass filtering. E.g. 2000"

  minCells:
    type: int?
    inputBinding:
      position: 7
      prefix: --minCells
    default: 1000
    doc: "The minimum number of cells in the ArchR project that must pass filtering before a warning message is printed. E.g. 1000"

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

  peaks_csv:
    type: File
    outputBinding:
      glob: "peaks.csv"

  peaks_bed:
    type: File
    outputBinding:
      glob: "peaks.bed"

  umap_coords_clusters_csv:
    type: File
    outputBinding:
      glob: "umap_coords_clusters.csv"

baseCommand: [Rscript, /opt/run_ArchR_analysis_pt2.R]