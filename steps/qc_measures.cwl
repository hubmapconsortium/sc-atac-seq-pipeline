#!/usr/bin/env cwl-runner

class: CommandLineTool
id: qc_measures
label: alignment/etc. QC
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.1.3
  InitialWorkDirRequirement:
    listing:
      - $(inputs.bam_file)
      - $(inputs.bam_index)


inputs:
  bam_file:
    type: File
    inputBinding:
      position: 0
      valueFrom: $(self.basename)

  bam_index:
    type: File
    inputBinding:
      position: 1
      valueFrom: $(self.basename)

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
