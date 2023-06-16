#!/usr/bin/env cwl-runner

class: CommandLineTool
id: sort_bam_file_tool
label: sort bam file
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.0.5

inputs:
  unsorted_paired_end_bam:
    type: File
    inputBinding:
      position: 0

  sorted_bam_name:
    type: string
    inputBinding:
      position: 1
    default: "sorted_alignment.bam"

  threads:
    type: int?
    inputBinding:
      position: 2
      prefix: "--threads"
    default: 0

outputs:
  sorted_paired_end_bam:
    type: File
    outputBinding:
      glob: $(inputs.sorted_bam_name)

baseCommand: /opt/bulk/sort_bam_file.py
