#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_preprocess_reads
label: snaptools preprocess reads
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

inputs:
  input_file:
    type: File
    inputBinding:
      position: 1
      prefix: --input-file
    doc: The genome BAM or BED file to be processed.

  output_snap:
    type: string?
    inputBinding:
      position: 2
      prefix: --output-snap
    default: "snaptools.snap"
    doc: The name to use for the output SNAP file.

  genome_name:
    type: string?
    inputBinding:
      position: 3
      prefix: --genome-name
    default: "hg38"
    doc: The genome identifier (i.e. hg19, hg38, mm10).

  genome_size:
    type: File
    inputBinding:
      position: 4
      prefix: --genome-size
    # TODO: why doesn't this default javascript work?
    #default: "$(inputs.input_file.basename).chrom.sizes"
    doc: A a text file that contains corresponding genome sizes.

  barcode_file:
    type: string?
    inputBinding:
      prefix: --barcode-file
      position: 6
    doc: A text file containing pre-selected cell barcodes.

  min_mapq:
    type: int?
    inputBinding:
      position: 7
      prefix: --min-mapq
    default: 30
    doc: Minimum mappability score.

  min_flen:
    type: int?
    inputBinding:
      position: 8
      prefix: --min-flen
    default: 0
    doc: Minimum fragment length.

  max_flen:
    type: int?
    inputBinding:
      position: 9
      prefix: --max-flen
    default: 1000
    doc: Maximum fragment length.

  min_cov:
    type: int?
    inputBinding:
      position: 10
      prefix: --min-cov
    default: 100
    doc: Minimum number of fragments per barcode.

  keep_chrm:
    type: string?
    inputBinding:
      position: 11
      prefix: --keep-chrm
    default: "TRUE"
    doc: A boolen tag indicating whether to keep fragments mapped to chrM.

  keep_single:
    type: string?
    inputBinding:
      position: 12
      prefix: --keep-single
    default: "FALSE"
    doc: A boolen tag indicating whether to keep those reads whose mates are not mapped or are missing.

  keep_secondary:
    type: string?
    inputBinding:
      position: 13
      prefix: --keep-secondary
    default: "FALSE"
    doc: A boolen tag indicating whether to keep secondary alignments.

  keep_discordant:
    type: string?
    inputBinding:
      position: 14
      prefix: --keep-discordant
    doc: A boolen tag indicating whether to keep discordant read pairs.

  tmp_folder:
    type: string?
    inputBinding:
      position: 15
      prefix: --tmp-folder
    default: "/tmp"
    doc: Directory to store temporary files.

  overwrite:
    type: string?
    inputBinding:
      position: 16
      prefix: --overwrite
    default: "TRUE"
    doc: Whether to overwrite the matrix session if it already exists.

  qc_file:
    type: string?
    inputBinding:
      position: 17
      prefix: --qc-file
    doc: A boolen tag indicating whether to create a master qc file.

  verbose:
    type: string?
    inputBinding:
      position: 18
      prefix: --verbose
    default: "TRUE"
    doc: A boolen tag; if true output the progress.

outputs:
  snap_file:
    type: File
    outputBinding:
      glob: $(inputs.output_snap)

  snap_qc_file:
    type: File
    outputBinding:
      glob: "$(inputs.output_snap).qc"


baseCommand: [snaptools, snap-pre]
