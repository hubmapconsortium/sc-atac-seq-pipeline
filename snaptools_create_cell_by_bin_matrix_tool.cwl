#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_create_cell_by_bin_matrix
label: snaptools create cell by bin matrix
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
      - entry: $(inputs.snap_file)
        writable: true

inputs:
  snap_file:
    type: File
    inputBinding:
      position: 1
      prefix: --snap-file
    doc: The SNAP file in which to place the cell by bin matrix.

  bin_size_list:
    type: int[]?
    inputBinding:
      prefix: --bin-size-list
      position: 2
    default: [5000, 10000]
    doc: A list of bin size(s) to create in the cell-by-bin count matrix.

  tmp_folder:
    type: string?
    inputBinding:
      position: 11
      prefix: --tmp-folder
    default: "/tmp"
    doc: Directory to store temporary files.

  verbose:
    type: string?
    inputBinding:
      position: 13
      prefix: --verbose
    default: "TRUE"
    doc: A boolen tag; if true output the progress.

outputs:
  snap_file_w_cell_by_bin:
    type: File
    outputBinding:
      glob: $(inputs.snap_file.basename)

baseCommand: [snaptools, snap-add-bmat]
