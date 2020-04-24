cwlVersion: v1.1
class: Workflow

# label: A workflow that processes bulk ATAC seq data
# doc: A workflow that processes bulk ATAC seq data

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
  if_sort: string?

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
      if_sort: if_sort

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
