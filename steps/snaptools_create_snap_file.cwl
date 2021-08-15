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

  fragment_file:
    type: File
    outputSource: snaptools_create_fragment_file/fragment_file

  analysis_CSV_files:
    type: File[]
    outputSource: analyze_with_ArchR/CSV_files

  analysis_PDF_files:
    type: File[]
    outputSource: analyze_with_ArchR/PDF_files

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


  align_reads:
    run: BWA-Mem.cwl
    in:
      alignment_index: snaptools_index_ref_genome/genome_alignment_index
      InputFile:
       source: adjust_barcodes/adj_fastq_dir
       valueFrom: |
          ${
            return [{"class":"File", "location": self.location + "/barcode_added_R1.fastq"},
                    {"class":"File", "location": self.location + "/barcode_added_R2.fastq"}]
          }

      # Index is provided by 'alignement_index input above
    out: [reads_stdout] 


  add_cell_identifiers_and_sort:
    run: add_cell_identifiers_and_sort.cwl
    in:
       sam_file: align_reads/reads_stdout
    out: [sorted_BAM_with_cell_ids]

  analyze_with_ArchR:
    run: analyze_snap_steps/archr_analyze.cwl
    in:
      bam_file: add_cell_identifiers_and_sort/sorted_BAM_with_cell_ids
      threads: threads
    out:
      - CSV_files
      - PDF_files

  snaptools_create_fragment_file:
    run: create_snap_steps/snaptools_create_fragment_file.cwl
    in:
      input_bam: add_cell_identifiers_and_sort/sorted_BAM_with_cell_ids
    out: [fragment_file]

