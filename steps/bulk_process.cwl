cwlVersion: v1.1
class: Workflow

# label: A workflow that processes bulk ATAC seq data
# doc: A workflow that processes bulk ATAC seq data

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  merged_bam: File
  alignment_index: Directory?
  encode_blacklist: File?
  threads: int?

outputs:
  bam_file:
    type: File
    outputSource: snaptools_remove_blacklist/rmsk_bam

  alignment_qc_report:
    type: File
    outputSource: alignment_qc/alignment_qc_report

steps:
  alignment_qc:
    run: create_snap_steps/alignment_qc.cwl
    in:
      bam_file: merged_bam
      threads: threads
    out: [alignment_qc_report]

  sort_bam_file:
    run: create_snap_steps/sort_bam_file_tool.cwl
    in:
      unsorted_paired_end_bam: merged_bam
      threads: threads
    out: [sorted_paired_end_bam]

  snaptools_remove_blacklist:
    run: create_snap_steps/snaptools_remove_blacklist.cwl
    in:
      bam_file: sort_bam_file/sorted_paired_end_bam
      bed_file: encode_blacklist
    out: [rmsk_bam]
