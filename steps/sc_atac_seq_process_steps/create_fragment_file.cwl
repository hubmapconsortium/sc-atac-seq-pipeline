#!/usr/bin/env cwl-runner

class: CommandLineTool
label: create fragment file
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.1.7

baseCommand: /opt/bulk/sort_index_frag.py

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
