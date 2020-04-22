cwlVersion: v1.1
class: Workflow

# label: A workflow that creates a SNAP file as outlined at: https://github.com/r3fang/SnapTools
# doc: A workflow that creates a SNAP file as outlined at: https://github.com/r3fang/SnapTools

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  input_reference_genome: File
  reference_genome_index: File
  genome_name: string?
  input_fastq1: File
  input_fastq2: File
  blacklist_bed: File?
  tmp_folder: string?
  alignment_threads: int?
  num_cores: int?

outputs:

  bam_file:
    type: File
    outputSource: snaptools_remove_blacklist/rmsk_bam

  zipped_files:
    type:
      type: array
      items: File
    outputSource: snaptools_fastqc_tool/zipped_files

  report_files:
    type:
      type: array
      items: File
    outputSource: snaptools_fastqc_tool/report_files

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

  motifs:
    type: File
    outputSource: motif_analysis/motifs

  known_results:
    type: File
    outputSource: motif_analysis/known_results

  autonormalization_statistics:
    type: File
    outputSource: motif_analysis/autonormalization_statistics

  formatted_results:
    type: File
    outputSource: motif_analysis/formatted_results


steps:

  snaptools_fastqc_tool:
    run: create_snap_steps/snaptools_fastqc_tool.cwl
    in:
      sequence_files: [input_fastq1, input_fastq2]
    out: [zipped_files, report_files]

  snaptools_align_paired_end:
    run: create_snap_steps/snaptools_align_paired_end_tool.cwl
    in:
      input_reference: reference_genome_index
      input_fastq1: input_fastq1
      input_fastq2: input_fastq2
      tmp_folder: tmp_folder
      num_threads: alignment_threads

    out: [paired_end_bam]

  sort_bam_file:
    run: create_snap_steps/sort_bam_file_tool.cwl
    in:
      unsorted_paired_end_bam: snaptools_align_paired_end/paired_end_bam
    out: [sorted_paired_end_bam]

  snaptools_remove_blacklist:
    run: create_snap_steps/snaptools_remove_blacklist.cwl
    in:
      bam_file: sort_bam_file/sorted_paired_end_bam
      bed_file: blacklist_bed
    out: [rmsk_bam]

  call_peaks:
    run: call_peaks.cwl
    in:
      bam_file: snaptools_remove_blacklist/rmsk_bam
    out: [peaks_table, narrow_peaks, summits_bed, r_script, bed_graphs]

  motif_analysis:
    run: motif_analysis.cwl
    in:
      bed_file: call_peaks/summits_bed
      num_cores: num_cores
    out: [motifs, known_results, autonormalization_statistics, formatted_results]
