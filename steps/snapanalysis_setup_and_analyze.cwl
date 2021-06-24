cwlVersion: v1.1
class: Workflow

inputs:
  encode_blacklist: File?
  gene_track: File?
  gene_annotation: File?
  preferred_barcodes: File?
  promoters: File?
  input_snap: File
  processes: int?

outputs:
  analysis_CSV_files:
    type:
      type: array
      items: File
    outputSource: snapanalysis_analyze/CSV_files

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

  peaks_bed_file:
    type: File
    outputSource: snapanalysis_analyze/peaks_bed_file

  umap_coords_csv:
    type: File
    outputSource: snapanalysis_analyze/umap_coords_csv

  motif_CSV_files:
    type:
      type: array
      items: File
    outputSource: snapanalysis_motif/CSV_files

  motif_RData_file:
    type: File
    outputSource: snapanalysis_motif/RData_file

  cell_by_bin_h5ad:
    type: File
    outputSource: convert_to_h5ad/cell_by_bin_h5ad

  cell_by_gene_h5ad:
    type: File
    outputSource: convert_to_h5ad/cell_by_gene_h5ad

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
      processes: processes
    out:
      - snap_rds
      - peaks_combined_bed
      - CSV_files
      - PDF_files
      - RDS_objects
      - peaks_bed_file
      - umap_coords_csv
      - cell_by_gene_raw_mtx
      - cell_by_gene_smoothed_hdf5
      - cell_by_bin_mtx
      - cell_by_bin_barcodes
      - cell_by_bin_bins

  snapanalysis_add_pmat_tool:
    run: analyze_snap_steps/snapanalysis_add_pmat_tool.cwl
    in:
      snap_file: input_snap
      peak_file: snapanalysis_analyze/peaks_combined_bed
    out:
      [snap_file_w_peaks]

  snapanalysis_motif:
    run: analyze_snap_steps/snapanalysis_motif.cwl
    in:
      snap_file: snapanalysis_add_pmat_tool/snap_file_w_peaks
      snap_rds: snapanalysis_analyze/snap_rds
    out:
      [CSV_files, RData_file]

  convert_to_h5ad:
    run: convert_to_h5ad.cwl
    in:
      umap_coords_csv: snapanalysis_analyze/umap_coords_csv
      cell_by_gene_raw_mtx: snapanalysis_analyze/cell_by_gene_raw_mtx
      cell_by_gene_smoothed_hdf5: snapanalysis_analyze/cell_by_gene_smoothed_hdf5
      cell_by_bin_mtx: snapanalysis_analyze/cell_by_bin_mtx
      cell_by_bin_barcodes: snapanalysis_analyze/cell_by_bin_barcodes
      cell_by_bin_bins: snapanalysis_analyze/cell_by_bin_bins
    out:
      - cell_by_bin_h5ad
      - cell_by_gene_h5ad
