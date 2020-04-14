cwlVersion: v1.1
class: Workflow

# label: A workflow that analyzes a SNAP file as outlined at: https://github.com/r3fang/SnapATAC
# doc: A workflow that analyzes a SNAP file as outlined at: https://github.com/r3fang/SnapATAC

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0001-5173-4627
    s:email: jshands@ucsc.edu
    s:name: Walter Shands

s:codeRepository: https://github.com/wshands/SnapATAC/tree/feature/snap-analysis
s:dateCreated: "2019-12-03"
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
  encode_blacklist: File
  gene_track: File
  gene_annotation: File?
  preferred_barcodes: File?
  promoters: File?
  input_snap: File

outputs:
  analysis_CSV_files:
    type:
      type: array
      items: File
    outputSource: snapanalysis_analyze/CSV_files

  analysis_BED_files:
    type:
      type: array
      items: File
    outputSource: snapanalysis_analyze/BED_files

  analysis_PDF_files:
    type:
      type: array
      items: File
    outputSource: snapanalysis_analyze/PDF_files

  analysis_RDS_objects:
    type:
      type: array
      items: File
    outputSource: snapanalysis_analyze/RDS_objects

  analysis_MTX_files:
    type:
      type: array
      items: File
    outputSource: snapanalysis_analyze/MTX_files

#  analysis_motif_file:
#    type: File
#    outputSource: snapanalysis_motif/motif_file


steps:
  snapanalysis_select_barcode:
    run: analyze_snap_steps/snapanalysis_select_barcode.cwl
    in:
      preferred_barcodes: preferred_barcodes
    out:
      [selected_barcodes, barcode_PDF_files]

  snapanalysis_analyze:
    run: analyze_snap_steps/snapanalysis_analyze.cwl
    in:
      input_snap: input_snap
      selected_barcodes: snapanalysis_select_barcode/selected_barcodes
      encode_blacklist: encode_blacklist
      gene_track: gene_track
      gene_annotation: gene_annotation
      promoters: promoters
    out:
      [peaks_combined_bed, CSV_files, BED_files, PDF_files, RDS_objects, MTX_files]

#  snapanalysis_add_pmat_tool:
#    run: analyze_snap_steps/snapanalysis_add_pmat_tool.cwl
#    in:
#      snap_file: input_snap
#      peak_file: snapanalysis_analyze/peaks_combined_bed
#    out:
#      [snap_file_w_peaks]
#
#  snapanalysis_motif:
#    run: analyze_snap_steps/snapanalysis_motif.cwl
#    in:
#      snap_file: snapanalysis_add_pmat_tool/snap_file_w_peaks
#    out:
#      [motif_file]
