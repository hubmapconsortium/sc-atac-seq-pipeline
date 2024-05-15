cwlVersion: v1.2
class: Workflow

# label: A workflow that processes and analyzes single cell ATAC seq data
# doc: A workflow that processes and analyzes single cell ATAC seq data

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  reference_genome_fasta: File?
  alignment_index: File?
  size_index: File?
  genome_name: string?
  sequence_directory: Directory[]
  tmp_folder: string?
  encode_blacklist: File?
  threads: int?
  if_sort: string?

outputs:

  fastqc_dir:
    type: Directory[]
    outputSource: fastqc/fastqc_dir

  bam_file:
    type: File
    outputSource: bulk_process/sorted_bam_file

  peaks_table:
    type: File
    outputSource: bulk_analysis/peaks_table

  narrow_peaks:
    type: File
    outputSource: bulk_analysis/narrow_peaks

  summits_bed:
    type: File
    outputSource: bulk_analysis/summits_bed

  r_script:
    type: File
    outputSource: bulk_analysis/r_script

  bed_graphs:
    type:
      type: array
      items: File
    outputSource: bulk_analysis/bed_graphs

  qc_report:
    type: File
    outputSource: qc_measures/qc_report

steps:
  fastqc:
    scatter: [fastq_dir]
    scatterMethod: dotproduct
    run: steps/fastqc.cwl
    in:
      fastq_dir: sequence_directory
      threads: threads
    out:
      [fastqc_dir]

  gather_sequence_bundles:
    run: bulk_gather_sequence_bundles.cwl
    in:
      sequence_directory: sequence_directory
    out:
      [fastq1_files, fastq2_files]

  align_paired_end:
    scatter: [input_fastq1, input_fastq2]
#    scatter: [Fastq_1, Fastq_2]
    scatterMethod: dotproduct
    run: steps/align_reads.cwl
    in:
#      alignment_index: index_ref_genome/genome_alignment_index
#      InputFile: [gather_sequence_bundles/fastq1_files, gather_sequence_bundles/fastq2_files]
#        source: [gather_sequence_bundles/fastq1_files, gather_sequence_bundles/fastq2_files]
#        linkMerge: merge_nested

      input_fastq1: gather_sequence_bundles/fastq1_files
      input_fastq2: gather_sequence_bundles/fastq2_files
#      tmp_folder: tmp_folder
#      num_threads: threads
#      if_sort: if_sort

    out: [paired_end_bam]

  merge_bam:
    run: steps/merge_bam.cwl
    in:
      bam_files: align_paired_end/paired_end_bam
    out: [merged_bam]

  index_merged_bam:
    run: steps/index_merged_bam.cwl
    in:
      merged_bam_file: merge_bam/merged_bam
    out:
      [sorted_merged_bam, merged_bam_index]

  bulk_process:
    run: steps/bulk_process.cwl
    in:
      sorted_merged_bam: index_merged_bam/sorted_merged_bam
      encode_blacklist: encode_blacklist
      threads: threads

    out:
      [sorted_bam_file]


  bulk_analysis:
    run: steps/bulk_analysis.cwl
    in:
      bam_file: bulk_process/sorted_bam_file

    out:
      [peaks_table, narrow_peaks, summits_bed, bed_graphs, r_script]

  qc_measures:
    run: steps/qc_measures.cwl
    in:
      bam_file: index_merged_bam/sorted_merged_bam
      peak_file: bulk_analysis/summits_bed
      bam_index: index_merged_bam/merged_bam_index
    out: [qc_report]
