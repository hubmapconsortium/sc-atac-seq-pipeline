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
    outputSource: remove_blacklist/rmsk_bam

steps:

  sort_bam_file:
    run: sc_atac_seq_process_steps/sort_bam_file_tool.cwl
    in:
      unsorted_paired_end_bam: merged_bam
      threads: threads
    out: [sorted_paired_end_bam]

  remove_blacklist:
    run: sc_atac_seq_process_steps/remove_blacklist.cwl
    in:
      bam_file: sort_bam_file/sorted_paired_end_bam
      bed_file: encode_blacklist
    out: [rmsk_bam]
