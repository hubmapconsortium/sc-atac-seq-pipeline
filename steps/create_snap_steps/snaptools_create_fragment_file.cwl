#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_create_fragment_file
label: snaptools create fragment file
cwlVersion: v1.1

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2020-03-30"

s:keywords: edam:topic_0091 , edam:topic_0622

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
    dockerPull: "quay.io/wshands/sc-atac-seq:latest"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:

  input_bam:
    type: File
    inputBinding:
      position: 1
      prefix: --input-bam
      default: "snaptools_alignment.bam"
      doc: The bam file to use

  output_bed:
    type: string
    inputBinding:
      position: 2
      prefix: --output-bed
    default: "fragment_file.bed"
    doc: The name to use for the output fragment file

  aligner:
    type: string?
    inputBinding:
      position: 4
      prefix: --aligner
    default: "bwa"
    doc: The name of the aligner, e.g. 'bwa'.

  path_to_aligner:
    type: string?
    inputBinding:
      position: 5
      prefix: --path-to-aligner
    default: "/tools"
    doc: The file system path to the aligner.

  aligner_options:
    type: string[]?
    inputBinding:
      prefix: --aligner-options
      position: 6
    doc: List of strings indicating options you would like passed to aligner.

  read_fastq_command:
    type: string?
    inputBinding:
      position: 7
      prefix: --read-fastq-command
    doc: Command line to execute for each of the input files.

  min_cov:
    type: string?
    inputBinding:
      position: 8
      prefix: --min-cov
    doc: Minimum number of fragments per barcode.

  num_threads:
    type: string?
    inputBinding:
      position: 9
      prefix: --num-threads
    doc: The number of threads to use.

  if_sort:
    type: string?
    inputBinding:
      position: 10
      prefix: --if-sort
    doc: Whether to sort the bam file based on the read name.

  tmp_folder:
    type: string?
    inputBinding:
      position: 11
      prefix: --tmp-folder
    default: "/tmp"
    doc: Directory to store temporary files.

  overwrite:
    type: string?
    inputBinding:
      position: 12
      prefix: --overwrite
    doc: Whether to overwrite the output file if it already exists.

  verbose:
    type: string?
    inputBinding:
      position: 13
      prefix: --verbose
    doc: A boolean tag; if true output the progress.

outputs:

  fragment_file:
    type: File
    outputBinding:
      glob: $(inputs.output_bed)



baseCommand: [snaptools, create-fragment-file]
