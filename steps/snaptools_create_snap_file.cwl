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


inputs:
  input_reference_genome: File
  reference_genome_index: File?
  genome_name: string?
  input_fastq1: File
  input_fastq2: File
  input_barcode_fastq: File?
  blacklist_bed: File
  tmp_folder: string?

outputs:
  snap_file:
    type: File
    outputSource: snaptools_create_cell_by_bin_matrix/snap_file_w_cell_by_bin
  snap_qc_file:
    type: File
    outputSource: snaptools_preprocess_reads/snap_qc_file

steps:
  snaptools_index_ref_genome:
    run: create_snap_steps/snaptools_index_ref_genome_tool.cwl
    in:
      input_fasta: input_reference_genome
      reference_genome_index: reference_genome_index

    out:
      [output]

  snaptools_create_ref_genome_size_file:
    run: create_snap_steps/snaptools_create_ref_genome_size_file_tool.cwl
    in:
      ref_genome: input_reference_genome

    out: [output]

  snaptools_add_barcodes_to_reads_tool:
    run: create_snap_steps/snaptools_add_barcodes_to_reads_tool.cwl
    in:
      input_fastq1: input_fastq1
      input_fastq2: input_fastq2
      input_barcode_fastq: input_barcode_fastq
     
    out: [barcode_added_fastq1, barcode_added_fastq2]

  snaptools_align_paired_end:
    run: create_snap_steps/snaptools_align_paired_end_tool.cwl
    in:
      input_reference: snaptools_index_ref_genome/output
      input_fastq1: snaptools_add_barcodes_to_reads_tool/barcode_added_fastq1
      input_fastq2: snaptools_add_barcodes_to_reads_tool/barcode_added_fastq2
      tmp_folder: tmp_folder

    out: [output]

  snaptools_remove_blacklist:
    run: create_snap_steps/snaptools_remove_blacklist.cwl
    in:
      bam_file: snaptools_align_paired_end/output
      bed_file: blacklist_bed

    out: [rmsk_bam]

  snaptools_preprocess_reads:
    run: create_snap_steps/snaptools_preprocess_reads_tool.cwl
    in:
      input_file: snaptools_remove_blacklist/rmsk_bam
      genome_size: snaptools_create_ref_genome_size_file/output
      genome_name: genome_name
      
    out: [snap_file, snap_qc_file]

  snaptools_create_cell_by_bin_matrix:
    run: create_snap_steps/snaptools_create_cell_by_bin_matrix_tool.cwl
    in:
      snap_file: snaptools_preprocess_reads/snap_file

    out: [snap_file_w_cell_by_bin]


