#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_index_ref_genome
label: snaptools index reference genome
cwlVersion: v1.1

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2019-11-15"
s:license: https://spdx.org/licenses/Apache-2.0

s:keywords: edam:topic_0091 , edam:topic_0622
s:programmingLanguage: Python

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

$schemas:
  - https://schema.org/docs/schema_org_rdfa.html
  - http://edamontology.org/EDAM_1.18.owl


dct:creator:
  '@id':  https://orcid.org/0000-0001-5173-4627
  foaf:name: Walter Shands
  foaf:mbox: jshands@ucsc.edu

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-grch38:1.0.1
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

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
