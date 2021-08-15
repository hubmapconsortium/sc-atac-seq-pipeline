#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_analyze
label: snap analysis analyze
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38
  NetworkAccess:
    networkAccess: true

inputs:
  bam-file:
    type: File
    inputBinding:
      position: 1
      prefix: --bam_file
    doc: The selected barcodes object file.

  processes:
    type: int?
    inputBinding:
      position: 7
      prefix: --processes
    default: 1
    doc: Number of processes to use

outputs:
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

  umap_coords_csv:
    type: File
    outputBinding:
      glob: "archr_umap_coords_clusters.csv"

baseCommand: [Rscript, /opt/run_ArchR_analysis.R]
