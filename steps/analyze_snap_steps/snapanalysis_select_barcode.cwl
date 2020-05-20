#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_select_barcode
label: snapanalysis select barcode
cwlVersion: v1.0

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
    dockerPull: "hubmap/sc-atac-seq-grch38"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:
  preferred_barcodes:
    type: File?
    inputBinding:
      position: 2
      prefix: --preferred_barcodes
    doc: A CSV file.

outputs:
  # An RDS output will not be produced
  # if no barcode CVS is provided as input
  selected_barcodes:
    type: ["null", File]
    outputBinding:
      glob: "*.rds"

  barcode_PDF_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.pdf"


baseCommand: [Rscript, /opt/snapAnalysis_select_barcode.R]
