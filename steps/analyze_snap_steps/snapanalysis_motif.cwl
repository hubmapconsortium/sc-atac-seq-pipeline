#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_motif
label: snap analysis motif
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.3
  NetworkAccess:
    networkAccess: true

inputs:
  snap_rds:
    type: File
    inputBinding:
      position: 1
      prefix: --snap_rds
    doc: The snap object RDS file.

  snap_file:
    type: File
    inputBinding:
      position: 2
      prefix: --snap_file
    doc: The snap file.

  #tmpdir:
  #  type: string?
  #  inputBinding:
  #    position: 7
  #    prefix: --tmpdir
  #  doc: Path for the temporary directory


outputs:
  CSV_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.csv"

  RData_file:
    type: File
    outputBinding:
      glob: "chromvar_data.RData"

baseCommand: [Rscript, /opt/snapMotifAnalysis.R]
