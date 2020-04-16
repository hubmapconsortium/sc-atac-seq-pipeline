cwlVersion: v1.1
class: Workflow

# label: A workflow that creates a SNAP file as outlined at: https://github.com/r3fang/SnapTools
# doc: A workflow that creates a SNAP file as outlined at: https://github.com/r3fang/SnapTools

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  input_reference_genome: File
  reference_genome_index: File?
  genome_name: string?
  input_fastq1: File
  input_fastq2: File
  input_barcode_fastq: File?
  blacklist_bed: File?
  tmp_folder: string?
  alignment_threads: int?
  if_sort: string?

outputs:
  bam_file:
    type: File
    outputSource: snaptools_remove_blacklist/rmsk_bam
  snap_file:
    type: File
    outputSource: snaptools_create_cell_by_bin_matrix/snap_file_w_cell_by_bin
  snap_qc_file:
    type: File
    outputSource: snaptools_preprocess_reads/snap_qc_file
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
  snaptools_index_ref_genome:
    run: create_snap_steps/snaptools_index_ref_genome_tool.cwl
    in:
      input_fasta: input_reference_genome
      reference_genome_index: reference_genome_index
    out:
      [genome_index]

  snaptools_create_ref_genome_size_file:
    run: create_snap_steps/snaptools_create_ref_genome_size_file_tool.cwl
    in:
      ref_genome: input_reference_genome
    out: [genome_sizes]

  snaptools_fastqc_tool:
    run: create_snap_steps/snaptools_fastqc_tool.cwl
    in:
      sequence_files: [input_fastq1, input_fastq2]
    out: [zipped_files, report_files]

  snaptools_align_paired_end:
    run: create_snap_steps/snaptools_align_paired_end_tool.cwl
    in:
      input_reference: snaptools_index_ref_genome/genome_index
      input_fastq1: input_fastq1
      input_fastq2: input_fastq2
      tmp_folder: tmp_folder
      num_threads: alignment_threads
      if_sort: if_sort

    out: [paired_end_bam]

  snaptools_remove_blacklist:
    run: create_snap_steps/snaptools_remove_blacklist.cwl
    in:
      bam_file: snaptools_align_paired_end/paired_end_bam
      bed_file: blacklist_bed
    out: [rmsk_bam]

  snaptools_preprocess_reads:
    run: create_snap_steps/snaptools_preprocess_reads_tool.cwl
    in:
      input_file: snaptools_remove_blacklist/rmsk_bam
      genome_size: snaptools_create_ref_genome_size_file/genome_sizes
      genome_name: genome_name
    out: [snap_file, snap_qc_file]

  snaptools_create_cell_by_bin_matrix:
    run: create_snap_steps/snaptools_create_cell_by_bin_matrix_tool.cwl
    in:
      snap_file: snaptools_preprocess_reads/snap_file
    out: [snap_file_w_cell_by_bin]


  #  sort_bam_file:
  #    run: create_snap_steps/sort_bam_file_tool.cwl
  #    in:
  #      unsorted_paired_end_bam: snaptools_align_paired_end/paired_end_bam
  #    out: [sorted_paired_end_bam]
