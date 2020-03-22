#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gather_sequence_bundles
label: gather sequence bundles
cwlVersion: v1.1

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2020-03-20"
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
    dockerPull: "quay.io/wshands/sc-atac-seq:feature_add-motif-analysis"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000
  InlineJavascriptRequirement: {}

  InitialWorkDirRequirement:
    listing:
      - $(inputs.sequence_directory)

inputs:
  sequence_directory:
    type: Directory
    inputBinding:
      position: 1
      #prefix: --sequence-directory
    doc: The directory with sample fastq or fastq.gz files.

outputs:
#  sequence_bundles:
#    type: File[]
#    outputBinding:
#      glob: "input.json"
#      loadContents: true
#      outputEval: |
#                ${
#                    var new_sequence_str = self[0].contents;
#                    var sequence_array_of_bundles = JSON.parse(new_sequence_str);
#                    var fastq2_array = [];
#                    for (var i =0; i < sequence_array_of_bundles.length; i++) {
#                        var fastq2 = sequence_array_of_bundles[0].input_fastq2;
#                        fastq2_array.push(fastq2);
#                    }
#                    return fastq2_array;
#                    //var fastq2_array_str = JSON.stringify(fastq2_array);
#                    //return fastq2_array_str;
#                 }

  fastq1_files:
    type: File[]
    outputBinding:
      glob: "input.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i =0; i < bundle_array.length; i++) {
                        var file = bundle_array[i].input_fastq1;
                        file_array.push(file);
                    }
                    return file_array;
                 }

  fastq2_files:
    type: File[]
    outputBinding:
      glob: "input.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i =0; i < bundle_array.length; i++) {
                        var file = bundle_array[i].input_fastq2;
                        file_array.push(file);
                    }
                    return file_array;
                 }

  barcode_fastq_files:
    type: File[]
    outputBinding:
      glob: "input.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i =0; i < bundle_array.length; i++) {
                        var file = bundle_array[i].input_barcode_fastq;
                        file_array.push(file);
                    }
                    return file_array;
                 }


baseCommand: [gather_sequence_files.py]
