cwlVersion: v1.1
class: Workflow

# label: A workflow that analyzes processed bulk ATAC seq data
# doc: A workflow that analyzes processed bulk ATAC seq data

inputs:
  bam_file:
    type: File

outputs:

  peaks_table:
    type: File
    outputSource: call_peaks/peaks_table

  narrow_peaks:
    type: File
    outputSource: call_peaks/narrow_peaks

  summits_bed:
    type: File
    outputSource: call_peaks/summits_bed

  r_script:
    type: File
    outputSource: call_peaks/r_script

  bed_graphs:
    type:
      type: array
      items: File
    outputSource: call_peaks/bed_graphs

steps:

  call_peaks:
    run: call_peaks.cwl
    in:
      bam_file: bam_file
    out: [peaks_table, narrow_peaks, summits_bed, r_script, bed_graphs]
