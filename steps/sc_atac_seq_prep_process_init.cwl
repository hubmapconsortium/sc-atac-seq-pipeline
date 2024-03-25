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
    type: File
    outputSource: sc_atac_seq_initial_analysis/bam_file

  bam_index:
    type: File
    outputSource: sc_atac_seq_initial_analysis/bam_index

  gene_row_data_csv:
    type: File
    outputSource: sc_atac_seq_initial_analysis/gene_row_data_csv

  cell_column_data_csv:
    type: File
    outputSource: sc_atac_seq_initial_analysis/cell_column_data_csv

  cell_by_bin_h5ad:
    type: File
    outputSource: sc_atac_seq_initial_analysis/cell_by_bin_h5ad

  cell_by_gene_h5ad:
    type: File
    outputSource: sc_atac_seq_initial_analysis/cell_by_gene_h5ad

  fragment_file:
    type: File
    outputSource: sc_atac_seq_initial_analysis/fragment_file

  image_file:
    type: File
    outputSource: sc_atac_seq_initial_analysis/image_file
  
  archr_project:
    type: Directory
    outputSource: sc_atac_seq_initial_analysis/archr_project

steps:
  fastqc:
    scatter: [fastq_dir]
    scatterMethod: dotproduct
    run: fastqc.cwl
    in:
      fastq_dir: sequence_directory
      threads: threads
    out:
      [fastqc_dir]

  concat_fastq:
    run: concat-fastq.cwl
    in:
      sequence_directory: sequence_directory
      assay: assay
    out:
      [output_directory, merged_fastq_r1, merged_fastq_r2, merged_fastq_barcode]

  sc_atac_seq_initial_analysis:
    run: sc_atac_seq_initial_analysis.cwl
    in:
      assay: assay
      concat_fastq_dir: concat_fastq/output_directory

      input_fastq1: concat_fastq/merged_fastq_r1
      input_fastq2: concat_fastq/merged_fastq_r2

      threads: threads
    out:
      - bam_file
      - bam_index
      - gene_row_data_csv
      - cell_column_data_csv
      - image_file
      - archr_project
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad
      - fragment_file