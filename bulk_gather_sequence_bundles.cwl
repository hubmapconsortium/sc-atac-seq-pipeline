#!/usr/bin/env cwl-runner

class: CommandLineTool
id: gather_sequence_bundles
label: gather sequence bundles
cwlVersion: v1.1

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

$schemas:
  - https://schema.org/docs/schema_org_rdfa.html
  - http://edamontology.org/EDAM_1.18.owl

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.1
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
    type: Directory[]
    inputBinding:
      position: 1
      #prefix: --sequence-directory
    doc: The directory with sample fastq or fastq.gz files.

outputs:
  fastq1_files:
    type: File[]
    outputBinding:
      # This file is written by gather_sequence_files.py
      glob: "input.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i =0; i < bundle_array.length; i++) {
                        // The dictionary key 'input_fastq1' must match the key used
                        // in the string template SEQUENCES_TEMPLATE in gather_sequence_files.py
                        var file = bundle_array[i].fastq_r1;
                        file_array.push(file);
                    }
                    return file_array;
                 }

  fastq2_files:
    type: File[]
    outputBinding:
      # This file is written by gather_sequence_files.py
      glob: "input.json"
      loadContents: true
      outputEval: |
                ${
                    var bundle_array_str = self[0].contents;
                    var bundle_array = JSON.parse(bundle_array_str);
                    var file_array = [];
                    for (var i =0; i < bundle_array.length; i++) {
                        // The dictionary key 'input_fastq2' must match the key used
                        // in the string template SEQUENCES_TEMPLATE in gather_sequence_files.py
                        var file = bundle_array[i].fastq_r2;
                        file_array.push(file);
                    }
                    return file_array;
                 }

baseCommand: [/opt/bulk/bulk_gather_sequence_files.py]
