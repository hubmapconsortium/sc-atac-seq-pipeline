cwlVersion: v1.2
class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  bam_file: File
  bam_index: File
  image_file: File
  threads: int?

outputs:

  Fragment_Size_Distribution_pdf:
    type: File
    outputSource: analyze_with_ArchR/Fragment_Size_Distribution_pdf

  TSS_by_Unique_Frags_pdf:
    type: File
    outputSource: analyze_with_ArchR/TSS_by_Unique_Frags_pdf

  QC-Sample-FragSizes-TSSProfile_pdf:
    type: File
    outputSource: analyze_with_ArchR/QC-Sample-FragSizes-TSSProfile_pdf

  QC-Sample-Statistics_pdf:
    type: File
    outputSource: analyze_with_ArchR/QC-Sample-Statistics_pdf

  TSS-vs-Frags_pdf:
    type: File
    outputSource: analyze_with_ArchR/TSS-vs-Frags_pdf

  Peak-Call-Summary_pdf:
    type: File
    outputSource: analyze_with_ArchR/Peak-Call-Summary_pdf

  Plot-UMAP-Sample-Clusters_pdf:
    type: File
    outputSource: analyze_with_ArchR/Plot-UMAP-Sample-Clusters_pdf

  peaks_bed:
    type: File
    outputSource: analyze_with_ArchR/peaks_bed

  umap_coords_clusters_csv:
    type: File
    outputSource: analyze_with_ArchR/umap_coords_clusters_csv

steps:

  analyze_with_ArchR:
    run: sc_atac_seq_analyze_steps/archr_clustering.cwl
    in:
      image_file
    out:
      - Fragment_Size_Distribution_pdf
      - TSS_by_Unique_Frags_pdf
      - QC-Sample-FragSizes-TSSProfile_pdf
      - QC-Sample-Statistics_pdf
      - TSS-vs-Frags_pdf
      - Peak-Call-Summary_pdf
      - Plot-UMAP-Sample-Clusters_pdf
      - peaks_bed
      - umap_coords_clusters_csv

baseCommand: [Rscript, /opt/run_ArchR_analysis_pt2.R]
