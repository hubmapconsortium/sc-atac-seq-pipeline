#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snapanalysis_analyze
label: snap analysis analyze
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
    dockerPull: hubmap/sc-atac-seq-grch38:1.1.2-sciseq
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000
  NetworkAccess:
    networkAccess: true
  # Set the environment variable for the location of
  # the R environment file. In this file will be set
  # the location of the R temporary working directory
  # to a directory on the host and not in
  # the Docker container. If this is not done R
  # can run out of memory due to memory limitations
  # in the container
  # https://stackoverflow.com/questions/27302983/how-to-change-directory-for-temporary-files-problems-with-huge-temporary-raste
  # https://stat.ethz.ch/R-manual/R-devel/library/base/html/Startup.html
  # https://www.commonwl.org/user_guide/12-env/index.html
  #EnvVarRequirement:
  #  envDef:
  #    #R_ENVIRON  : $(inputs.r_environ_path)
  #    TEMP  : $(inputs.tmpdir)
  #    TMP  : $(inputs.tmpdir)
  #    TMPDIR  : $(inputs.tmpdir)

inputs:
  input_snap:
    type: File
    inputBinding:
      position: 1
      prefix: --input_snap
    doc: The selected barcodes object file.

  selected_barcodes:
    type: File?
    inputBinding:
      position: 2
      prefix: --selected_barcodes
    doc: The selected barcodes object file.

  encode_blacklist:
    type: File?
    inputBinding:
      position: 3
      prefix: --encode_blacklist
    doc: A BED file of ENCODE blacklist to prevent potential artifacts.

  gene_track:
    type: File?
    inputBinding:
      position: 4
      prefix: --gene_track
    doc: A BED file of gene tracks.

  gene_annotation:
    type: File?
    inputBinding:
      position: 5
      prefix: --gene_annotation
    doc: A GTF file of gene annotations.

  promoters:
    type: File?
    inputBinding:
      position: 6
      prefix: --promoters
    doc: A BED file of promoters.

  processes:
    type: int?
    inputBinding:
      position: 7
      prefix: --processes
    default: 1
    doc: Number of processes to use

  #tmpdir:
  #  type: string?
  #  inputBinding:
  #    position: 7
  #    prefix: --tmpdir
  #  doc: Path for the temporary directory


outputs:
  snap_rds:
    type: File
    outputBinding:
      glob: "peaks_snap.rds"

  peaks_combined_bed:
    type: File
    outputBinding:
      glob: "peaks.combined.bed"

  CSV_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.csv"

  BED_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.bed"

  PDF_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.pdf"

  RDS_objects:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.rds"

  TXT_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.txt"

  MTX_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.mtx"

  HDF5_files:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.hdf5"

baseCommand: [Rscript, /opt/snapAnalysis.R]
