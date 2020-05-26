cwlVersion: v1.1
class: Workflow

# label: A workflow that processes and analyzes single cell ATAC seq data
# doc: A workflow that processes and analyzes single cell ATAC seq data

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  reference_genome_fasta: File?
  alignment_index: File?
  size_index: File?
  genome_name: string?
  sequence_directory: Directory
  blacklist_bed: File?
  tmp_folder: string?

  encode_blacklist: File
  threads: int?
  if_sort: string?

outputs:

  fastqc_dir:
    type: Directory
    outputSource: fastqc/fastqc_dir

  bam_file:
    type:
      type: array
      items: File
    outputSource: bulk_process/bam_file

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

steps:

  fastqc:
    run: fastqc.cwl
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

  bulk_process:
    scatter: [input_fastq1, input_fastq2]
    scatterMethod: dotproduct
    run: steps/bulk_process.cwl
    in:
      reference_genome_fasta: reference_genome_fasta
      alignment_index: alignment_index
      size_index: size_index
      genome_name: genome_name
      input_fastq1: gather_sequence_bundles/fastq1_files
      input_fastq2: gather_sequence_bundles/fastq2_files
      blacklist_bed: blacklist_bed
      tmp_folder: tmp_folder
      threads: threads
      if_sort: if_sort

    out:
      [bam_file]


  bulk_analysis:
    run: steps/bulk_analysis.cwl
    in:
      bam_files: bulk_process/bam_file

    out:
      [peaks_table, narrow_peaks, summits_bed, bed_graphs, r_script]
