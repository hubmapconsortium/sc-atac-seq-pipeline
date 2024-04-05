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
    outputSource: sc_atac_seq_prep_process_init/fastqc_dir

  bam_file:
    type: File?
    outputSource: maybe_save_bam_file/bam_output

  fragment_file:
    type: File
    outputSource: sc_atac_seq_prep_process_init/fragment_file

  Fragment_Size_Distribution_pdf:
    type: File
    outputSource: sc_atac_seq_prep_process_init/Fragment_Size_Distribution_pdf

  TSS_by_Unique_Frags_pdf:
    type: File
    outputSource: sc_atac_seq_prep_process_init/TSS_by_Unique_Frags_pdf

  QC-Sample-FragSizes-TSSProfile_pdf:
    type: File
    outputSource: sc_atac_seq_prep_process_init/QC-Sample-FragSizes-TSSProfile_pdf

  QC-Sample-Statistics_pdf:
    type: File
    outputSource: sc_atac_seq_prep_process_init/QC-Sample-Statistics_pdf

  TSS-vs-Frags_pdf:
    type: File
    outputSource: sc_atac_seq_prep_process_init/TSS-vs-Frags_pdf

  Peak-Call-Summary_pdf:
    type: File
    outputSource: analyze_with_ArchR/Peak-Call-Summary_pdf

  Plot-UMAP-Sample-Clusters_pdf:
    type: File
    outputSource: analyze_with_ArchR/Plot-UMAP-Sample-Clusters_pdf

  gene_row_data_csv:
    type: File
    outputSource: sc_atac_seq_prep_process_init/gene_row_data_csv

  cell_column_data_csv:
    type: File
    outputSource: sc_atac_seq_prep_process_init/cell_column_data_csv

  umap_coords_clusters_csv:
    type: File
    outputSource: analyze_with_ArchR/umap_coords_clusters_csv


  cell_by_bin_h5ad:
    type: File
    outputSource: sc_atac_seq_prep_process_init/cell_by_bin_h5ad

  cell_by_gene_h5ad:
    type: File
    outputSource: sc_atac_seq_prep_process_init/cell_by_gene_h5ad

  genome_build_json:
    type: File
    outputSource: write_genome_build/genome_build_json

steps:

  sc_atac_seq_prep_process_init:
    run: steps/sc_atac_seq_prep_process_init.cwl
    in:
     assay: assay
     concat_fastq_dir: concat_fastq/output_directory
     orig_fastq_dir: sequence_directory
     input_fastq1: concat_fastq/merged_fastq_r1
     input_fastq2: concat_fastq/merged_fastq_r2

     threads: threads
    out:
      - fastqc_dir
      - bam_file
      - bam_index
      - gene_row_data_csv
      - cell_column_data_csv
      - fragment_file
      - TSS_by_Unique_Frags_pdf
      - Fragment_Size_Distribution_pdf
      - image_file
      - archr_project
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad
      - QC-Sample-FragSizes-TSSProfile_pdf
      - QC-Sample-Statistics_pdf
      - TSS-vs-Frags_pdf

  analyze_with_ArchR:
    run: steps/sc_atac_seq_analyze_steps/archr_clustering.cwl
    in:
      image_file: sc_atac_seq_prep_process_init/image_file
      archr_project: sc_atac_seq_prep_process_init/archr_project
    out:
      - Peak-Call-Summary_pdf
      - Plot-UMAP-Sample-Clusters_pdf
      - peaks_bed
      - umap_coords_clusters_csv

  qc_measures:
    run: steps/qc_measures.cwl
    in:
      bam_file: sc_atac_seq_prep_process_init/bam_file
      bam_index: sc_atac_seq_prep_process_init/bam_index
      peak_file: analyze_with_ArchR/peaks_bed
      cell_by_bin_h5ad: sc_atac_seq_prep_process_init/cell_by_bin_h5ad
    out:
      - qc_report


  write_genome_build:
    run: steps/write_genome_build.cwl
    in: {}
    out: [genome_build_json]

  #thanks to @pvanheus in the CWL gitter instance
  maybe_save_bam_file:
    in:
      bam_input: sc_atac_seq_prep_process_init/bam_file
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
