cwlVersion: v1.1
class: Workflow

# label: A workflow that creates a SNAP file as outlined at: https://github.com/r3fang/SnapTools
# doc: A workflow that creates a SNAP file as outlined at: https://github.com/r3fang/SnapTools

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapTools/tree/feature/docker_cwl
s:dateCreated: "2019-11-15"
s:license: https://spdx.org/licenses/Apache-2.0

s:keywords: edam:topic_0091 , edam:topic_0622
s:programmingLanguage: Python

$namespaces:
  s: https://schema.org/
  edam: http://edamontology.org/

$schemas:
  - https://schema.org/docs/schema_org_rdfa.html
  - http://edamontology.org/EDAM_1.18.owl

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
  blacklist_bed: File?
  tmp_folder: string?
  threads: int?
  processes: int?
  bin_size_list: int[]?

  encode_blacklist: File?
  gene_track: File?
  gene_annotation: File?
  preferred_barcodes: File?
  promoters: File?


outputs:
  bam_file:
    type: File
    outputSource: snaptools_remove_blacklist/rmsk_bam
  alignment_qc_report:
    type: File
    outputSource: alignment_qc/alignment_qc_report
  fragment_file:
    type: File
    outputSource: snaptools_create_fragment_file/fragment_file
  snap_file:
    type: File
    outputSource: snaptools_create_cell_by_bin_matrix/snap_file_w_cell_by_bin
  snap_qc_file:
    type: File
    outputSource: snaptools_preprocess_reads/snap_qc_file

  analysis_CSV_files:
    type: File[]
    outputSource: snapanalysis_setup_and_analyze/analysis_CSV_files

  analysis_BED_files:
    type: File[]
    outputSource: snapanalysis_setup_and_analyze/analysis_BED_files

  analysis_PDF_files:
    type: File[]
    outputSource: snapanalysis_setup_and_analyze/analysis_PDF_files

  analysis_RDS_objects:
    type: File[]
    outputSource: snapanalysis_setup_and_analyze/analysis_RDS_objects

  cell_by_gene_matrix:
    type: File
    outputSource: snapanalysis_setup_and_analyze/cell_by_gene_matrix

  cell_by_bin_mtx:
    type: File
    outputSource: snapanalysis_setup_and_analyze/cell_by_bin_mtx

  cell_by_bin_barcodes:
    type: File
    outputSource: snapanalysis_setup_and_analyze/cell_by_bin_barcodes

  cell_by_bin_bins:
    type: File
    outputSource: snapanalysis_setup_and_analyze/cell_by_bin_bins

  motif_CSV_files:
    type: File[]
    outputSource: snapanalysis_setup_and_analyze/motif_CSV_files

  motif_RData_file:
    type: File
    outputSource: snapanalysis_setup_and_analyze/motif_RData_file

steps:
  snaptools_index_ref_genome:
    run: create_snap_steps/snaptools_index_ref_genome_tool.cwl
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

  snaptools_align_paired_end:
    run: create_snap_steps/snaptools_align_paired_end_tool.cwl
    in:
      alignment_index: snaptools_index_ref_genome/genome_alignment_index

      input_fastq1:
        source: adjust_barcodes/adj_fastq_dir
        valueFrom: |
           ${
             return {"class":"File", "location": self.location + "/barcode_added_R1.fastq"}
           }

      input_fastq2:
        source: adjust_barcodes/adj_fastq_dir
        valueFrom: |
           ${
             return {"class":"File", "location": self.location + "/barcode_added_R2.fastq"}
           }

      tmp_folder: tmp_folder
      num_threads: threads

    out: [paired_end_bam]

  alignment_qc:
    run: create_snap_steps/alignment_qc.cwl
    in:
      bam_file: snaptools_align_paired_end/paired_end_bam
      threads: threads
    out: [alignment_qc_report]

  snaptools_remove_blacklist:
    run: create_snap_steps/snaptools_remove_blacklist.cwl
    in:
      bam_file: snaptools_align_paired_end/paired_end_bam
      bed_file: blacklist_bed
    out: [rmsk_bam]

  snaptools_create_fragment_file:
    run: create_snap_steps/snaptools_create_fragment_file.cwl
    in:
      input_bam: snaptools_align_paired_end/paired_end_bam
    out: [fragment_file]

  snaptools_preprocess_reads:
    run: create_snap_steps/snaptools_preprocess_reads_tool.cwl
    in:
      input_bam: snaptools_remove_blacklist/rmsk_bam
      genome_size: snaptools_index_ref_genome/genome_size_index
      genome_name: genome_name
    out: [snap_file, snap_qc_file]

  snaptools_create_cell_by_bin_matrix:
    run: create_snap_steps/snaptools_create_cell_by_bin_matrix_tool.cwl
    in:
      snap_file: snaptools_preprocess_reads/snap_file
      bin_size_list: bin_size_list
    out: [snap_file_w_cell_by_bin]

  snapanalysis_setup_and_analyze:
    run: snapanalysis_setup_and_analyze.cwl
    in:
      input_snap: snaptools_create_cell_by_bin_matrix/snap_file_w_cell_by_bin
      preferred_barcodes: preferred_barcodes
      encode_blacklist: encode_blacklist
      gene_track: gene_track
      gene_annotation: gene_annotation
      promoters: promoters

    out:
      - analysis_CSV_files
      - analysis_BED_files
      - analysis_PDF_files
      - analysis_RDS_objects
      - cell_by_gene_matrix
      - cell_by_bin_mtx
      - cell_by_bin_barcodes
      - cell_by_bin_bins
      - motif_CSV_files
      - motif_RData_file
