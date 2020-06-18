#!/usr/bin/env cwl-runner

class: CommandLineTool
id: alignment_qc_tool
label: alignment qc
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.0
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:
  bam_file:
    type: File
    inputBinding:
      position: 0

  threads:
    type: int?
    inputBinding:
      position: 1
    default: 1

outputs:
  alignment_qc_report:
    type: File
    outputBinding:
      glob: alignment_qc.json

baseCommand: /opt/alignment_qc.py
