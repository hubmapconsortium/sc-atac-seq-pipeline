cwlVersion: v1.2
class: CommandLineTool
label: Assay-specific adjustment of cell barcodes
hints:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-barcode-adj:2.1.4
baseCommand: /opt/adjust_barcodes.py

inputs:
  assay:
    type: string
    inputBinding:
      position: 0
  directory:
    type: Directory[]
    inputBinding:
      prefix: --directory
      position: 1
  orig_dir:
    type: Directory[]
    inputBinding:
      prefix: --orig_dir
  metadata_file:
    type: File?
    inputBinding:
      prefix: --metadata_file
outputs:
  adj_fastq_dir:
    type: Directory
    outputBinding:
      glob: 'adj_fastq'
  metadata_json:
    type: File?
    outputBinding:
      glob: 'metadata.json'

