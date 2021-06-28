#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_add_barcodes_to_reads_tool
label: snaptools add barcodes to reads
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.1
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

inputs:
  input_fastq1:
    type: File
    inputBinding:
      position: 1
      prefix: --input-fastq1
    doc: The first paired end fastq file.

  input_fastq2:
    type: File
    inputBinding:
      position: 2
      prefix: --input-fastq2
    doc: The second paired end fastq file.

  input_barcode_fastq:
    type: File?
    inputBinding:
      position: 3
      prefix: --input-barcode-fastq
    doc: The fastq file containing the cell barcodes.

  output_fastq_prefix:
    type: string?
    inputBinding:
      position: 4
      prefix: --output-fastq-prefix
    default: "barcode_added"
    doc: The prefix to use for the output fastq files with barcodes added to the read sequence identifier.


outputs:
  barcode_added_fastq1:
    type: File
    outputBinding:
      glob: $(inputs.output_fastq_prefix).R1.fastq

  barcode_added_fastq2:
    type: File
    outputBinding:
      glob: $(inputs.output_fastq_prefix).R3.fastq

baseCommand: [/opt/add_barcodes_to_reads.pl]
