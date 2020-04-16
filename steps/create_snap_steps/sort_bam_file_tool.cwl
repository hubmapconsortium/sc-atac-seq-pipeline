#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_align_paired_end
label: snaptools align paired end reads
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: "seandonahue5311/bulk-atac:v1.0"
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

baseCommand: /opt/sort_bam_file.py
