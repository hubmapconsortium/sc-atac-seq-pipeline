cwlVersion: v1.2
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
hints:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-barcode-adj:2.0b1
baseCommand: /opt/adjust_barcodes.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  directory:
    type: Directory[]
    inputBinding:
      position: 1
outputs:
  adj_fastq_dir:
    type: Directory
    outputBinding:
      glob: 'adj_fastq'
  metadata_json:
    type: File?
    outputBinding:
      glob: 'metadata.json'

