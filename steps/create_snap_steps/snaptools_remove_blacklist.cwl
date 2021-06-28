#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_remove_blacklist
label: snaptools remove blacklist
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.4.1

inputs:
  bam_file:
    type: File
    inputBinding:
      position: 1
      prefix: --bam-file
    doc: The genome BAM file to be processed.

  bed_file:
    type: File?
    inputBinding:
      position: 2
      prefix: --bed-file
    doc: The genome BED file with the blacklisted regions to be removed.


outputs:
  rmsk_bam:
    type: File
    outputBinding:
      glob: rmsk.bam

baseCommand: [/opt/remove_blacklist.sh]
