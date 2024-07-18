#!/usr/bin/env cwl-runner
class: CommandLineTool
id: archr_clustering
label: ArchR clustering
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.1.5
    dockerOutputDirectory: "/output"
  NetworkAccess:
    networkAccess: true
  InitialWorkDirRequirement:
    listing:
      - $(inputs.archr_project)

inputs:
  image_file:
    type: File
    inputBinding:
      position: 1
      prefix: --image_file
    doc: "The R image from the previous archR analysis step"

  archr_project:
    type: Directory
    inputBinding:
      position: 2
      prefix: --archr_project
      valueFrom: $(self.basename)
    doc: "The ArchRProj directory from the previous step"

outputs:
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
