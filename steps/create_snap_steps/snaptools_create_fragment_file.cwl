#!/usr/bin/env cwl-runner

class: CommandLineTool
label: snaptools create fragment file
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: "seandonahue5311/atac:v0.2"

baseCommand: /opt/sort_index_frag.py

inputs:
  input_bam:
    type: File
    inputBinding:
      position: 1
    doc: The bam file to use

outputs:
  fragment_file:
    type: File
    outputBinding:
      glob: "*.sort.bed.gz"
