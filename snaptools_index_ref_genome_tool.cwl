#!/usr/bin/env cwl-runner

class: CommandLineTool
id: snaptools_index_ref_genome
label: snaptools index reference genome
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
      - $(inputs.input_fasta)

inputs:
  reference_genome_index:
    type: File?
    inputBinding:
      position: 1
      prefix: --reference-genome-index
    doc: The genome reference index file in tar.gz format. If provided snaptools index-genome is not called.

  input_fasta:
    type: File
    inputBinding:
      position: 2
      prefix: --input-fasta
    doc: The genome reference file to be indexed.

  output_prefix:
    type: string?
    inputBinding:
      position: 3
      prefix: --output-prefix
    #default: "$(inputs.input_fasta[\"nameroot\"])"
    doc: The output prefix for the genome reference index file.

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

  num_threads:
    type: string?
    inputBinding:
      position: 6
      prefix: --num-threads
    doc: The number of threads to use.

outputs:
  output:
    type: File
    secondaryFiles: [".bwt", ".sa", ".ann", ".pac", ".amb"]
    outputBinding:
      glob: $(inputs.input_fasta.basename)

baseCommand: [create_reference_genome_index.sh]
