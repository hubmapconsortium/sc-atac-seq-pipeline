#!/usr/bin/env cwl-runner

class: CommandLineTool
id: call_peaks.cwl
label: call peaks
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.0.5
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:

  bam_file:
    type: File
    inputBinding:
      position: 1
      prefix: -t
    doc: The bam file to call peaks on

outputs:

  peaks_table:
    type: File
    outputBinding:
      glob: "*.xls"
    doc: Tabular file containing peak data

  narrow_peaks:
    type: File
    outputBinding:
      glob: "*.narrowPeak"
    doc: BED6+4 format file which contains the peak locations together with peak summit, p-value, and q-value. Can be loaded into UCSC genome browser

  summits_bed:
    type: File
    outputBinding:
      glob: "*.bed"
    doc: BED file containing summit locations for every peak

  r_script:
    type: File
    outputBinding:
      glob: "*.r"
    doc: An r script for generating a pdf of model based on data

  bed_graphs:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.bdg"
    doc: Bed graph files containing information about pileup signals and local biases


baseCommand: [macs2, callpeak]
