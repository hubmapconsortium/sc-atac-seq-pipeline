#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_add_pmat_tool
label: snapanalysis add pmat
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
    dockerPull: "quay.io/wshands/sc-atac-seq:feature_initial-pipeline"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.snap_file)
        writable: true


inputs:
  snap_file:
    type: File
    inputBinding:
      position: 1
      prefix: --snap-file
    doc: The SNAP file to be processed.

  peak_file:
    type: File
    inputBinding:
      position: 2
      prefix: --peak-file
    doc: The BED file that contains peak information.

  buffer_size:
    type: string?
    inputBinding:
      position: 3
      prefix: --buffer-size
    doc: Max number of barcodes be stored in memory (default 1000).

  tmp_folder:
    type: string?
    inputBinding:
      position: 4
      prefix: --tmp-folder
    default: "/tmp"
    doc: Directory to store temporary files.

  verbose:
    type: string?
    inputBinding:
      position: 5
      prefix: --verbose
    default: "TRUE"
    doc: A boolen tag; if true output the progress.

outputs:
  snap_file_w_peaks:
    type: File
    outputBinding:
      glob: $(inputs.snap_file)

baseCommand: [snaptools, snap-add-pmat]
