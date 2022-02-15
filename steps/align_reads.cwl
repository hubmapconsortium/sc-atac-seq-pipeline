#!/usr/bin/env cwl-runner

class: CommandLineTool
id: align_reads
label: align paired end reads with HISAT2
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/hubmap/sc-atac-hisat2-hg38:latest

inputs:
  input_fastq1:
    type: File
    inputBinding:
      position: 3
    doc: The first paired end fastq file to be aligned.

  input_fastq2:
    type: File
    inputBinding:
      position: 4
    doc: The second paired end fastq file to be aligned.

  num_threads:
    type: int?
    inputBinding:
      position: 1
      prefix: --processes
    default: 16
    doc: The number of threads to use.

outputs:
  paired_end_bam:
    type: File
    outputBinding:
      glob: alignment.bam

baseCommand: [/opt/align_reads.py]
