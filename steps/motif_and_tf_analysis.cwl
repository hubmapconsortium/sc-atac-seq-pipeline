#This will be a command line tool which wraps the r script that Matt wrote
#!/usr/bin/env cwl-runner

class: CommandLineTool
id: motif_and_tf_analysis.cwl
label: motif and tf analysis
cwlVersion: v1.1

requirements:
  DockerRequirement:
    dockerPull: hubmap/sc-atac-seq-hg38:2.0.6
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
    outdirMin: 100000

inputs:

  narrow_peaks:
    type: File
    inputBinding:
      position: 1
    doc: narrow peaks file

  bam_files:
    type:
      type: array
      items: File
    inputBinding:
      position: 2
    doc: A list of bam files

outputs:

  motifs:
    type: File
    outputBinding:
      glob:

  scores_matrix:
    type: File
    outputBinding:
      glob:

baseCommand: [Rscript, /opt/chromvar-analysis.R]
