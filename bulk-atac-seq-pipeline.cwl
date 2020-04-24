cwlVersion: v1.1
class: Workflow

# label: A workflow that processes and analyzes single cell ATAC seq data
# doc: A workflow that processes and analyzes single cell ATAC seq data

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  input_reference_genome: File
  reference_genome_index: File?
  genome_name: string?
  sequence_directory: Directory
  blacklist_bed: File?
  tmp_folder: string?

  encode_blacklist: File
  alignment_threads: int?
  if_sort: string?

outputs:

  zipped_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_bulk_snap_file/zipped_files

  report_files:
    type:
      type: array
      items:
         type: array
         items: File
    outputSource: create_bulk_snap_file/report_files

  bam_file:
    type: File[]
    outputSource: create_bulk_snap_file/bam_file

  peaks_table:
    type: File
    outputBinding:
      glob: "*.xls"
    doc: Tabular file containing peak data

  narrow_peaks:
    type: File
    outputBinding:
      glob: "*.narrowPeak"
    doc: BED6+4 format file which contains the peak locations together with peak summit, p-value, and q-value. Can be loaded into UCSC genome browser

  summits_bed:
    type: File
    outputBinding:
      glob: "*.bed"
    doc: BED file containing summit locations for every peak

  r_script:
    type: File
    outputBinding:
      glob: "*.r"
    doc: An r script for generating a pdf of model based on data

  bed_graphs:
    type:
      type: array
      items: File
    outputBinding:
      glob: "*.bdg"
    doc: Bed graph files containing information about pileup signals and local biases

steps:
  gather_sequence_bundles:
    run: bulk_gather_sequence_bundles.cwl
    in:
      sequence_directory: sequence_directory
    out:
      [fastq1_files, fastq2_files]

  bulk_process:
    scatter: [input_fastq1, input_fastq2]
    scatterMethod: dotproduct
    run: steps/bulk_process.cwl
    in:
     input_reference_genome: input_reference_genome
     reference_genome_index: reference_genome_index
     genome_name: genome_name
     input_fastq1: gather_sequence_bundles/fastq1_files
     input_fastq2: gather_sequence_bundles/fastq2_files
     blacklist_bed: blacklist_bed
     tmp_folder: tmp_folder
     alignment_threads: alignment_threads
     if_sort: if_sort

    out:
      [zipped_files, report_files, bam_file]


  bulk_analysis:
    run: steps/bulk_analysis.cwl
    in:
      bam_files: bulk_process/bam_file

    out:
      #motifs
      #scores
