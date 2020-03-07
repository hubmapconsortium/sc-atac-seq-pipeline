#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_add_barcodes_to_reads_tool
label: snaptools add barcodes to reads
cwlVersion: v1.1

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2019-01-09"
s:license: https://spdx.org/licenses/Apache-2.0

s:keywords: edam:topic_0091 , edam:topic_0622
s:programmingLanguage: Pearl

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

inputs:
  input_fastq1:
    type: File
    inputBinding:
      position: 1
      prefix: --input-fastq1
    doc: The first paired end fastq file.

  input_fastq2:
    type: File
    inputBinding:
      position: 2
      prefix: --input-fastq2
    doc: The second paired end fastq file.

  input_barcode_fastq:
    type: File
    inputBinding:
      position: 3
      prefix: --input-barcode-fastq
    doc: The fastq file containing the cell barcodes.

  output_fastq_prefix:
    type: string?
    inputBinding:
      position: 4
      prefix: --output-fastq-prefix
    default: "barcode_added"
    doc: The prefix to use for the output fastq files with barcodes added to the read sequence identifier.

outputs:
  barcode_added_fastq1:
    type: File
    outputBinding:
      glob: $(inputs.output_fastq_prefix).R1.fastq

  barcode_added_fastq2:
    type: File
    outputBinding:
      glob: $(inputs.output_fastq_prefix).R3.fastq

baseCommand: [add_barcodes_to_reads.pl]
