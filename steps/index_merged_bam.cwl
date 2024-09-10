#!/usr/bin/env cwl-runner

class: CommandLineTool
id: align_reads
label: align paired end reads with HISAT2
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-hisat2-hg38:2.1.7

inputs:
  merged_bam_file:
    type: File
    inputBinding:
      position: 3
    doc: Merged BAM file

  num_threads:
    type: int?
    inputBinding:
      position: 1
      prefix: --processes
    default: 16
    doc: The number of threads to use.

outputs:
  sorted_merged_bam:
    type: File
    outputBinding:
      glob: alignment.bam
  merged_bam_index:
    type: File
    outputBinding:
      glob: alignment.bam.bai

baseCommand: [/opt/index_merged_bam.py]
