#!/usr/bin/env cwl-runner

class: CommandLineTool
id: motif_analysis.cwl
label: motif analysis
cwlVersion: v1.0

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/homer:4.9.1--pl5.22.0_1"
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:

  bed_file:
    type: File
    inputBinding:
      position: 0
      prefix: -t
    doc: The bed file to analyze motifs from

  genome_name:
    type: string?
    inputBinding:
      position: 1
    default: "hg38"
    doc: The name of the genome being used

  output_directory:
    type: string?
    inputBinding:
      position: 2
    default: "."
    doc: The directory for Homer to write output files to

  region_size:
    type: int?
    inputBinding:
      position: 3
      prefix: -size
    default: 200
    doc: Size of the region around peak summit for motif finding

  rna_flag:
    type: string?
    inputBinding:
      position: 4
    default: "-rna"
    doc: A flag to let Homer know to search for RNA motifs

  num_cores:
    type: int?
    inputBinding:
      position: 5
      prefix: -p
    default: 1
    doc: The number of cores to use for motif analysis


outputs:

  motifs:
    type: File
    outputBinding:
      glob: "*.all.motifs"
    doc: A file containing all Homer motifs found

  known_results:
    type: File
    outputBinding:
      glob: "*Results.txt"
    doc: text file containing statistics about known motif enrichment

  autonormalization_statistics:
    type: File
    outputBinding:
      glob: "*autonorm.tsv"
    doc: autonormalization statistics for lower-order oligo normalization.

  formatted_results:
    type: File
    outputBinding:
      glob: "*Results.html"
    doc: formatted output of de novo motif finding.

baseCommand: [/usr/local/share/homer-4.9.1-1/bin/findMotifsGenome.pl]
