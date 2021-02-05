#!/usr/bin/env cwl-runner
class: CommandLineTool
id: merge_bam.cwl
label: Merge alignment output
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.3

arguments:
  - merged.bam

inputs:
  bam_files:
    type: File[]
    inputBinding:
      position: 1

outputs:
  merged_bam:
    type: File
    outputBinding:
      glob: "merged.bam"

baseCommand: [samtools, merge]
