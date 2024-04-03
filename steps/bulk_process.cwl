cwlVersion: v1.1
class: Workflow

# label: A workflow that processes bulk ATAC seq data
# doc: A workflow that processes bulk ATAC seq data

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  sorted_merged_bam: File
  encode_blacklist: File?
  threads: int?

outputs:
  sorted_bam_file:
    type: File
    outputSource: remove_blacklist/rmsk_bam

steps:

  remove_blacklist:
    run: sc_atac_seq_process_steps/remove_blacklist.cwl
    in:
      bam_file: sorted_merged_bam
      bed_file: encode_blacklist
    out: [rmsk_bam]
