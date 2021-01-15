#!/usr/bin/env cwl-runner

class: CommandLineTool
id: concat_fastq
label: Concatenate FASTQ files
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38
  InlineJavascriptRequirement: {}

inputs:
  sequence_directory:
    type: Directory[]
    inputBinding:
      position: 1
      #prefix: --sequence-directory
    doc: A directory with sample fastq or fastq.gz files.

  assay:
    type: string
    inputBinding:
      position: 2
    doc: The assay type of the input files.

outputs:
  output_directory:
    type: Directory
    outputBinding:
      glob: "concat_output_dir"
 
  merged_fastq_r1:
    type: File
    outputBinding:
      glob: "concat_output_dir/merged_R1.fastq"

  merged_fastq_r2:
    type: File
    outputBinding:
      glob: "concat_output_dir/merged_R2.fastq"

  merged_fastq_barcode:
    type: File?
    outputBinding:
      glob: "concat_output_dir/merged_R3.fastq"
      outputEval: ${ if((self === null || self.length == 0)) { return null;} else { return self[0].basename="merged_barcode.fastq"; } }

baseCommand: [/opt/concat_fastq.py]
