#!/usr/bin/env cwl-runner

class: CommandLineTool
id: concat_fastq
label: Concatenate FASTQ files
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.1-snare
  InlineJavascriptRequirement: {}

inputs:
  sequence_directory:
    type: Directory[]
    inputBinding:
      position: 1
      #prefix: --sequence-directory
    doc: A directory with sample fastq or fastq.gz files.

outputs:
  merged_fastq_r1:
    type: File
    outputBinding:
      glob: "merged_R1.fastq"

  merged_fastq_r2:
    type: File
    outputBinding:
      glob: "merged_R2.fastq"

  merged_fastq_barcode:
    type: File
    outputBinding:
      glob: "merged_barcode.fastq"

baseCommand: [/opt/concat_fastq.py]
