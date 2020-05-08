cwlVersion: v1.1
class: Workflow

# label: A workflow that creates and analyzes a SNAP file as outlined at:
# https://github.com/r3fang/SnapTools and https://github.com/r3fang/SnapATAC
# doc: A workflow that analyzes a SNAP file as outlined at: https://github.com/r3fang/SnapATAC

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapATAC/tree/feature/snap-analysis
s:dateCreated: "2020-02-15"
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
  reference_genome_fasta: File?
  alignment_index: File?
  size_index: File?
  genome_name: string?
  sequence_directory: Directory
  blacklist_bed: File?
  tmp_folder: string?

  encode_blacklist: File?
  gene_track: File?
  gene_annotation: File?
  preferred_barcodes: File?
  promoters: File?

  alignment_threads: int?

outputs:
  zipped_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/zipped_files

  report_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/report_files

  bam_file:
    type: File[]
    outputSource: create_and_analyze_snap_file/bam_file

  snap_file:
    type: File[]
    outputSource: create_and_analyze_snap_file/snap_file

  snap_qc_file:
    type: File[]
    outputSource: create_and_analyze_snap_file/snap_qc_file

  analysis_motif_file:
    type: File[]
    outputSource: create_and_analyze_snap_file/analysis_motif_file

  analysis_CSV_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/analysis_CSV_files

  analysis_BED_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/analysis_BED_files

  analysis_PDF_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/analysis_PDF_files

  analysis_RDS_objects:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/analysis_RDS_objects

  analysis_TXT_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/analysis_TXT_files

  analysis_MTX_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_and_analyze_snap_file/analysis_MTX_files

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

steps:
  gather_sequence_bundles:
    run: gather_sequence_bundles.cwl
    in:
      sequence_directory: sequence_directory
    out:
      [fastq1_files, fastq2_files, barcode_fastq_files]

  create_and_analyze_snap_file:
    scatter: [input_fastq1, input_fastq2, input_barcode_fastq]
    scatterMethod: dotproduct
    run: steps/snaptools_create_snap_file.cwl
    in:
     reference_genome_fasta: reference_genome_fasta
     alignment_index: alignment_index
     size_index: size_index
     genome_name: genome_name
     input_fastq1: gather_sequence_bundles/fastq1_files
     input_fastq2: gather_sequence_bundles/fastq2_files
     input_barcode_fastq: gather_sequence_bundles/barcode_fastq_files
     blacklist_bed: blacklist_bed
     tmp_folder: tmp_folder
     alignment_threads: alignment_threads

     preferred_barcodes: preferred_barcodes
     encode_blacklist: encode_blacklist
     gene_track: gene_track
     gene_annotation: gene_annotation
     promoters: promoters

    out:
      [zipped_files, report_files, bam_file, snap_file, snap_qc_file,
      analysis_motif_file, analysis_CSV_files, analysis_BED_files, analysis_PDF_files,
      analysis_RDS_objects, analysis_TXT_files, analysis_MTX_files]
