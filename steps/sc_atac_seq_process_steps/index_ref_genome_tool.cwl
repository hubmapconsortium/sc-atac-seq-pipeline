#!/usr/bin/env cwl-runner

class: CommandLineTool
id: index_ref_genome
label: index reference genome
cwlVersion: v1.2

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.0b7

  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta)

inputs:
  input_fasta:
    type: File?
    inputBinding:
      position: 1
      prefix: --reference-genome
    doc: The reference genome file in FASTA format, if in

  alignment_index:
    type: File?
    inputBinding:
      position: 2
      prefix: --alignment-index
    doc: The alignment index file in tar.* format. If provided snaptools index-genome is not called.

  size_index:
    type: File?
    inputBinding:
      position: 3
      prefix: --size-index
    doc: The genome size index, produced by "samtools faidx".

outputs:
  genome_alignment_index:
    type: Directory?
    outputBinding:
      glob: "index"
  genome_size_index:
    type: File?
    outputBinding:
      glob: "*.fai"

baseCommand: [/opt/index_reference_genome.py]
