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

inputs:
  image_file:
    type: File
    inputBinding:
      position: 3
      prefix: --image
    doc: "The R image from the previous archR analysis step"

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