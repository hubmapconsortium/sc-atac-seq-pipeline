#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_fastqc_tool
label: snaptools fastqc tool
cwlVersion: v1.0

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2020-03-19"
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
    dockerPull: "biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000
  InitialWorkDirRequirement:
    listing:
      - $(inputs.sequence_files)

inputs:
  outdir:
    type: string?
    inputBinding:
      position: 1
      prefix: --outdir
    doc: Create all output files in the specified output directory.

  extract:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --extract
    doc: Do not uncompress the output file after creating it.

  sequence_files:
    type: File[]
    inputBinding:
      position: 3
    doc: set of sequence files for each of which a quality control report is produced

outputs:
  zipped_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.zip"
    doc: Individual graph files and additional data files containing the raw data from which plots were drawn.
      
  report_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.html"
    doc: HTML reports with embedded graphs


baseCommand: [fastqc]
