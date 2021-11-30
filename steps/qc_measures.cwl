#!/usr/bin/env cwl-runner

class: CommandLineTool
id: qc_measures
label: alignment/etc. QC
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.4

inputs:
  bam_file:
    type: File
    inputBinding:
      position: 0

  peak_file:
    type: File
    inputBinding:
      position: 1

  cell_by_bin_h5ad:
    type: File?
    inputBinding:
      position: 2

outputs:
  qc_report:
    type: File
    outputBinding:
      glob: qc_report.json

baseCommand: /opt/qc_measures.py
