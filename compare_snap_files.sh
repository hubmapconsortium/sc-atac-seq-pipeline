#!/bin/bash

#Rscript compare_snap_files.R --first_barcode_csv /mnt/huBMAP_snapATAC_testing/Dinhs_files/cellBarcodes.csv \
#       	--second_barcode_csv /mnt/huBMAP_snapATAC_testing/Dinhs_files/cellBarcodes.csv \
#       	--first_cell_by_bin_summary_csv /mnt/huBMAP_snapATAC_testing/Dinhs_files/cellByBin_summary.csv \
#	--second_cell_by_bin_summary_csv /mnt/huBMAP_snapATAC_testing/Dinhs_files/cellByBin_summary.csv \
#	--barcode_prefix "KC20_" 
#
Rscript compare_snap_files.R --first_barcode_csv /mnt/huBMAP_snapATAC_testing/walts_run_using_Dinhs_input/cellBarcodes.csv \
       	--second_barcode_csv /mnt/huBMAP_snapATAC_testing/Dinhs_files/cellBarcodes.csv \
       	--first_cell_by_bin_summary_csv /mnt/huBMAP_snapATAC_testing/walts_run_using_Dinhs_input/cellByBin_summary.csv \
	--second_cell_by_bin_summary_csv /mnt/huBMAP_snapATAC_testing/Dinhs_files/cellByBin_summary.csv \
	--barcode_prefix "KC20_" 
