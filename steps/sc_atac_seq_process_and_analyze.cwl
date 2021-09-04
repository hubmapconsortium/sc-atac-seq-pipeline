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

  reference_genome_fasta: File?
  alignment_index: File?
  size_index: File?
  genome_name: string?
  input_fastq1: File
  input_fastq2: File
  threads: int?

outputs:
  unsorted_reads:
    type: File
    outputSource: align_reads/reads_stdout

  bam_file:
    type: File
    outputSource: add_cell_identifiers_and_sort/sorted_BAM_with_cell_ids

#  fragment_file:
#    type: File
#    outputSource: create_fragment_file/fragment_file

  analysis_CSV_files:
    type: File[]
    outputSource: analyze_with_ArchR/CSV_files

  analysis_PDF_files:
    type: File[]
    outputSource: analyze_with_ArchR/PDF_files


#  cell_by_bin_h5ad:
#    type: File
#    outputSource: convert_to_h5ad/cell_by_bin_h5ad
#
#  cell_by_gene_h5ad:
#    type: File
#    outputSource: convert_to_h5ad/cell_by_gene_h5ad


steps:
  index_ref_genome:
    run: sc_atac_seq_process_steps/index_ref_genome_tool.cwl
    in:
      input_fasta: reference_genome_fasta
      alignment_index: alignment_index
      size_index: size_index
    out:
      [genome_alignment_index, genome_size_index]

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
    run: BWA-Mem.cwl
    in:
      alignment_index: index_ref_genome/genome_alignment_index
#      InputFile:
#       source: adjust_barcodes/adj_fastq_dir
#       valueFrom: |
#          ${
#            return [{"class":"File", "location": self.location + "/barcode_added_R1.fastq"},
#                    {"class":"File", "location": self.location + "/barcode_added_R2.fastq"}]
#          }


      Fastq_1: 
        source: adjust_barcodes/adj_fastq_dir
        valueFrom: |
          ${
            return {"class":"File", "location": self.location + "/barcode_added_R1.fastq"}
          }

      Fastq_2:
        source: adjust_barcodes/adj_fastq_dir
        valueFrom: |
          ${
            return {"class":"File", "location": self.location + "/barcode_added_R2.fastq"}
          }



      # Index is provided by 'alignment_index input above
    out: [reads_stdout] 


  add_cell_identifiers_and_sort:
    run: add_cell_identifiers_and_sort.cwl
    in:
       sam_file: align_reads/reads_stdout
    out: [sorted_BAM_with_cell_ids]

  analyze_with_ArchR:
    run: sc_atac_seq_analyze_steps/archr_analyze.cwl
    in:
      bam_file: add_cell_identifiers_and_sort/sorted_BAM_with_cell_ids
      threads: threads
    out:
      - CSV_files
      - PDF_files

#  create_fragment_file:
#    run: sc_atac_seq_process_steps/create_fragment_file.cwl
#    in:
#      input_bam: add_cell_identifiers_and_sort/sorted_BAM_with_cell_ids
#    out: [fragment_file]



#  convert_to_h5ad:
#    run: convert_to_h5ad.cwl
#    in:
#      umap_coords_csv: snapanalysis_analyze/umap_coords_csv
#      cell_by_gene_raw_mtx: snapanalysis_analyze/cell_by_gene_raw_mtx
#      cell_by_gene_smoothed_hdf5: snapanalysis_analyze/cell_by_gene_smoothed_hdf5
#      cell_by_bin_mtx: snapanalysis_analyze/cell_by_bin_mtx
#      cell_by_bin_barcodes: snapanalysis_analyze/cell_by_bin_barcodes
#      cell_by_bin_bins: snapanalysis_analyze/cell_by_bin_bins
#    out:
#      - cell_by_bin_h5ad
#      - cell_by_gene_h5ad
