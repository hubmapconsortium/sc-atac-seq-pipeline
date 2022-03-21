#!/usr/bin/env cwl-runner

class: CommandLineTool
id: write_genome_build
label: Write genome build
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.0b7

inputs: {}

outputs:
  genome_build_json:
    type: File
    outputBinding:
      glob: genome_build.json

baseCommand: /opt/write_genome_build.py
