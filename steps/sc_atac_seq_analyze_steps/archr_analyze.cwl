#!/usr/bin/env cwl-runner

class: CommandLineTool
id: archr_analyze
label: ArchR analysis
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38
  NetworkAccess:
    networkAccess: true

arguments: ["--QCDir", $(runtime.outdir)]

inputs:
  bam_file:
    type: File
    inputBinding:
      position: 1
      prefix: --bam_file
    doc: The sorted BAM file with cell ids in the CB tag.
 
  threads:
    type: int?
    inputBinding:
      position: 3
      prefix: --threads
    default: 2
    doc: Number of threads to use

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

#  umap_coords_csv:
#    type: File
#    outputBinding:
#      glob: "archr_umap_coords_clusters.csv"

baseCommand: [Rscript, /opt/run_ArchR_analysis.R]
