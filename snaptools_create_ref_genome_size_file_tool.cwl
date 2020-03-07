#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_create_ref_genome_sizes_file
label: snaptools create reference genome sizes file
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
    dockerPull: "snaptools:1.2.3"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

  InitialWorkDirRequirement:
    listing:
      - $(inputs.ref_genome)

inputs:
  ref_genome:
    type: File
    inputBinding:
      position: 1
    doc: The genome reference file.

  output_file:
    type: string?
    inputBinding:
      position: 2
    default: "reference_genome_sizes.txt"
    doc: The name of the genome size output file.

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_file)

baseCommand: [create_genome_size_file.sh]
