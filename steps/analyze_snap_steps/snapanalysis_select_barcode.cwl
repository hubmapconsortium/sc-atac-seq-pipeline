#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_select_barcode
label: snapanalysis select barcode
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.4
  NetworkAccess:
    networkAccess: true

inputs:
  preferred_barcodes:
    type: File?
    inputBinding:
      position: 2
      prefix: --preferred_barcodes
    doc: A CSV file.

outputs:
  # An RDS output will not be produced
  # if no barcode CVS is provided as input
  selected_barcodes:
    type: ["null", File]
    outputBinding:
      glob: "*.rds"

  barcode_PDF_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.pdf"


baseCommand: [Rscript, /opt/snapAnalysis_select_barcode.R]
