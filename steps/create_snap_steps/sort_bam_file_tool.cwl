#!/usr/bin/env cwl-runner

class: CommandLineTool
id: sort_bam_file_tool
label: sort bam file
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:
  unsorted_paired_end_bam:
    type: File
    inputBinding:
      position: 0

  sorted_bam_name:
    type: string
    inputBinding:
      position: 1
    default: "sorted-snaptools_alignment.bam"

outputs:
  sorted_paired_end_bam:
    type: File
    outputBinding:
      glob: $(inputs.sorted_bam_name)

baseCommand: /opt/bulk/sort_bam_file.py
