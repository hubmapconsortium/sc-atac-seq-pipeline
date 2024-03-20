
cwlVersion: v1.2
class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  assay: string
  concat_fastq_dir: Directory

  input_fastq1: File
  input_fastq2: File
  threads: int?

outputs:
  bam_file:
    type: File
    outputSource: align_reads/paired_end_bam

  bam_index:
    type: File
    outputSource: align_reads/paired_end_bam_index

  cell_column_data_csv:
    type: File
    outputSource: analyze_with_ArchR/cell_column_data_csv

  gene_row_data_csv:
    type: File
    outputSource: analyze_with_ArchR/gene_row_data_csv
  
  fragment_file:
    type: File
    outputSource: create_fragment_file/fragment_file

  image_file:
    type: File
    outputSource: analyze_with_ArchR/image_file

  archr_project:
    type: Directory[]
    outputSource: analyze_with_ArchR/archr_project

  cell_by_bin_h5ad:
    type: File
    outputSource: convert_to_h5ad/cell_by_bin_h5ad

  cell_by_gene_h5ad:
    type: File
    outputSource: convert_to_h5ad/cell_by_gene_h5ad

steps:
  adjust_barcodes:
    run: adjust-barcodes.cwl
    in:
     assay: assay
     directory:
       # https://www.commonwl.org/user_guide/misc/ Connect a solo value to an input that expects an array of that type
       source: [ concat_fastq_dir ]
       linkMerge: merge_nested
    out:
     [adj_fastq_dir]

  align_reads:
    run: align_reads.cwl
    in:
      num_threads: threads

      input_fastq1:
        source: adjust_barcodes/adj_fastq_dir
        valueFrom: |
          ${
            return {"class": "File", "location": self.location + "/barcode_added_R1.fastq"}
          }

      input_fastq2:
        source: adjust_barcodes/adj_fastq_dir
        valueFrom: |
          ${
            return {"class": "File", "location": self.location + "/barcode_added_R2.fastq"}
          }

    out: [paired_end_bam, paired_end_bam_index]

  create_fragment_file:
    run: sc_atac_seq_process_steps/create_fragment_file.cwl
    in:
      input_bam: align_reads/paired_end_bam
    out: [fragment_file]

  analyze_with_ArchR:
    run: sc_atac_seq_analyze_steps/archr_init_analyze.cwl
    in:
      bam_file: align_reads/paired_end_bam
      bam_index: align_reads/paired_end_bam_index
      threads: threads
    out:
      - cell_column_data_csv
      - gene_row_data_csv
      - cell_by_gene_raw_mtx
      - cell_by_gene_smoothed_hdf5
      - cell_by_bin_mtx
      - cell_by_bin_barcodes
      - cell_by_bin_bins
      - image_file
      - archr_project

  convert_to_h5ad:
    run: convert_to_h5ad.cwl
    in:
      cell_column_data: analyze_with_ArchR/cell_column_data_csv
      cell_by_gene_raw_mtx: analyze_with_ArchR/cell_by_gene_raw_mtx
      cell_by_gene_smoothed_hdf5: analyze_with_ArchR/cell_by_gene_smoothed_hdf5
      cell_by_bin_mtx: analyze_with_ArchR/cell_by_bin_mtx
      cell_by_bin_barcodes: analyze_with_ArchR/cell_by_bin_barcodes
      cell_by_bin_bins: analyze_with_ArchR/cell_by_bin_bins
    out:
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad

