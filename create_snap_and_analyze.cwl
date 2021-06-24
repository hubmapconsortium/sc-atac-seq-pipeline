cwlVersion: v1.1
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  assay: string
  reference_genome_fasta: File?
  alignment_index: File?
  size_index: File?
  genome_name: string?
  sequence_directory: Directory[]
  blacklist_bed: File?
  tmp_folder: string?
  bin_size_list: int[]?

  encode_blacklist: File?
  gene_track: File?
  gene_annotation: File?
  preferred_barcodes: File?
  promoters: File?

  threads: int?
  exclude_bam: boolean?

outputs:
  fastqc_dir:
    type: Directory[]
    outputSource: fastqc/fastqc_dir

  fragment_file:
    type: File
    outputSource: create_and_analyze_snap_file/fragment_file

  bam_file:
    type: File?
    outputSource: maybe_save_bam_file/bam_output

  qc_report:
    type: File
    outputSource: qc_measures/qc_report

  snap_file:
    type: File
    outputSource: create_and_analyze_snap_file/snap_file

  snap_qc_file:
    type: File
    outputSource: create_and_analyze_snap_file/snap_qc_file

  analysis_CSV_files:
    type: File[]
    outputSource: create_and_analyze_snap_file/analysis_CSV_files

  analysis_PDF_files:
    type: File[]
    outputSource: create_and_analyze_snap_file/analysis_PDF_files

  analysis_RDS_objects:
    type: File[]
    outputSource: create_and_analyze_snap_file/analysis_RDS_objects

  motif_CSV_files:
    type: File[]
    outputSource: create_and_analyze_snap_file/motif_CSV_files

  motif_RData_file:
    type: File
    outputSource: create_and_analyze_snap_file/motif_RData_file

  cell_by_bin_h5ad:
    type: File
    outputSource: create_and_analyze_snap_file/cell_by_bin_h5ad

  cell_by_gene_h5ad:
    type: File
    outputSource: create_and_analyze_snap_file/cell_by_gene_h5ad

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

  concat_fastq:
    run: steps/concat-fastq.cwl
    in:
      sequence_directory: sequence_directory
      assay: assay
    out:
      [output_directory, merged_fastq_r1, merged_fastq_r2, merged_fastq_barcode]

  create_and_analyze_snap_file:
    run: steps/snaptools_create_snap_file.cwl
    in:
     assay: assay
     concat_fastq_dir: concat_fastq/output_directory

     reference_genome_fasta: reference_genome_fasta
     alignment_index: alignment_index
     size_index: size_index
     genome_name: genome_name

     input_fastq1: concat_fastq/merged_fastq_r1
     input_fastq2: concat_fastq/merged_fastq_r2


     blacklist_bed: blacklist_bed
     tmp_folder: tmp_folder
     threads: threads
     processes: threads
     bin_size_list: bin_size_list

     preferred_barcodes: preferred_barcodes
     encode_blacklist: encode_blacklist
     gene_track: gene_track
     gene_annotation: gene_annotation
     promoters: promoters

    out:
      - bam_file
      - fragment_file
      - snap_file
      - snap_qc_file
      - analysis_CSV_files
      - analysis_PDF_files
      - analysis_RDS_objects
      - peaks_bed_file
      - motif_CSV_files
      - motif_RData_file
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad

  qc_measures:
    run: steps/qc_measures.cwl
    in:
      bam_file: create_and_analyze_snap_file/bam_file
      peak_file: create_and_analyze_snap_file/peaks_bed_file
      cell_by_bin_h5ad: create_and_analyze_snap_file/cell_by_bin_h5ad
    out:
      - qc_report

  # thanks to @pvanheus in the CWL gitter instance
  maybe_save_bam_file:
    in:
      bam_input: create_and_analyze_snap_file/bam_file
      exclude_bam: exclude_bam
    out:
      - bam_output
    run:
      class: ExpressionTool
      inputs:
        bam_input:
          type: File
        exclude_bam:
          type: boolean?
          default: false
      outputs:
        bam_output:
          type: File?
      expression: |-
        ${
          if (inputs.exclude_bam) {
            return { 'bam_output': null };
          } else {
            return { 'bam_output': inputs.bam_input };
          }
        }
