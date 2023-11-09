cwlVersion: v1.2
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  assay: string
  sequence_directory: Directory[]

  threads: int?
  exclude_bam: boolean?

outputs:
  fastqc_dir:
    type: Directory[]
    outputSource: fastqc/fastqc_dir

  bam_file:
    type: File?
    outputSource: maybe_save_bam_file/bam_output

  fragment_file:
    type: File
    outputSource: sc_atac_seq_process_and_analyze/fragment_file

  Fragment_Size_Distribution_pdf:
    type: File
    outputSource: clustering/Fragment_Size_Distribution_pdf

  TSS_by_Unique_Frags_pdf:
    type: File
    outputSource: clustering/TSS_by_Unique_Frags_pdf

  QC-Sample-FragSizes-TSSProfile_pdf:
    type: File
    outputSource: clustering/QC-Sample-FragSizes-TSSProfile_pdf

  QC-Sample-Statistics_pdf:
    type: File
    outputSource: clustering/QC-Sample-Statistics_pdf

  TSS-vs-Frags_pdf:
    type: File
    outputSource: clustering/TSS-vs-Frags_pdf

  gene_row_data_csv:
    type: File
    outputSource: sc_atac_seq_process_and_analyze/gene_row_data_csv

  cell_column_data_csv:
    type: File
    outputSource: sc_atac_seq_process_and_analyze/cell_column_data_csv

  cell_by_bin_h5ad:
    type: File
    outputSource: sc_atac_seq_process_and_analyze/cell_by_bin_h5ad

  cell_by_gene_h5ad:
    type: File
    outputSource: sc_atac_seq_process_and_analyze/cell_by_gene_h5ad

  genome_build_json:
    type: File
    outputSource: write_genome_build/genome_build_json

steps:
#  fastqc:
#    scatter: [fastq_dir]
#    scatterMethod: dotproduct
#    run: steps/fastqc.cwl
#    in:
#      fastq_dir: sequence_directory
#      threads: threads
#    out:
#      [fastqc_dir]

#  concat_fastq:
#    run: steps/concat-fastq.cwl
#    in:
#      sequence_directory: sequence_directory
#      assay: assay
#    out:
#      [output_directory, merged_fastq_r1, merged_fastq_r2, merged_fastq_barcode]

  sc_atac_seq_process_and_analyze:
    run: steps/sc_atac_seq_process_and_analyze.cwl
    in:
      sequence_directory: sequence_directory
      assay: assay
      threads: threads
#      concat_fastq_dir: concat_fastq/output_directory

#      input_fastq1: concat_fastq/merged_fastq_r1
#      input_fastq2: concat_fastq/merged_fastq_r2

      threads: threads
    out:
      - bam_file
      - bam_index
      - fragment_file
      - gene_row_data_csv
      - cell_column_data_csv
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad
      - rdata_file

  clustering:
    run: steps/sc_atac_seq_analyze_steps/archr_clustering.cwl
    in:
      rdata_file: sc_atac_seq_process_and_analyze/rdata_file
    out:
      - Fragment_Size_Distribution_pdf
      - TSS_by_Unique_Frags_pdf
      - QC-Sample-FragSizes-TSSProfile_pdf
      - QC-Sample-Statistics_pdf
      - TSS-vs-Frags_pdf 

  qc_measures:
    run: steps/qc_measures.cwl
    in:
      bam_file: sc_atac_seq_process_and_analyze/bam_file
      bam_index: sc_atac_seq_process_and_analyze/bam_index
      cell_by_bin_h5ad: sc_atac_seq_process_and_analyze/cell_by_bin_h5ad
    out:
      - qc_report


  write_genome_build:
    run: steps/write_genome_build.cwl
    in: {}
    out: [genome_build_json]

  #thanks to @pvanheus in the CWL gitter instance
  maybe_save_bam_file:
    in:
      bam_input: sc_atac_seq_process_and_analyze/bam_file
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