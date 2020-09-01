#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_remove_blacklist
label: snaptools remove blacklist
cwlVersion: v1.1

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2020-3-4"
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
    dockerPull: hubmap/sc-atac-seq-grch38:1.1-snare
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

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
