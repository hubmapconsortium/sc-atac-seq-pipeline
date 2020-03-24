#!/usr/bin/env Rscript
install.packages("optparse")
library("optparse")

option_list = list(
  make_option(c("-b", "--first_barcode_csv"), type="character", default=NULL,
              help="First barcode file", metavar="character"),
  make_option(c("-c", "--second_barcode_csv"), type="character", default=NULL,
              help="Second barcode file", metavar="character"),
  make_option(c("-s", "--first_cell_by_bin_summary_csv"), type="character", default=NULL,
              help="First cell by bin summary csv", metavar="character"),
  make_option(c("-t", "--second_cell_by_bin_summary_csv"), type="character", default=NULL,
              help="Second cell by bin summary csv", metavar="character"),
  make_option(c("-p", "--barcode_prefix"), type="character", default=NULL,
              help="Prefix string on barcodes", metavar="character")

  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$first_barcode_csv)){
  print_help(opt_parser)
  stop("--first_barcode_csv argument must be supplied", call.=FALSE)
}
if (is.null(opt$second_barcode_csv)){
  print_help(opt_parser)
  stop("--second_barcode_csv argument must be supplied", call.=FALSE)
}


mytmpdir <- Sys.getenv("TMP")
sprintf("TMP is: %s", mytmpdir)
mytmpdir <- Sys.getenv("TMPDIR")
sprintf("TMPDIR is: %s", mytmpdir)
mytmpdir <- Sys.getenv("TEMP")
sprintf("TEMP is: %s", mytmpdir)
mytmpdir <- tempdir()
sprintf("tempdir() is: %s", mytmpdir)

message(sprintf("Reading first barcode csv\n"))
first_barcode_dataframe <- read.csv(opt$first_barcode_csv)
# Strip barcode prefix if it is there
first_barcode_dataframe$Barcode <- sub(paste0("^", opt$barcode_prefix), "", first_barcode_dataframe$Barcode)
head(first_barcode_dataframe)
nrow(first_barcode_dataframe)

message(sprintf("Reading second barcode csv\n"))
second_barcode_dataframe <- read.csv(opt$second_barcode_csv)
second_barcode_dataframe$Barcode <- sub(paste0("^", opt$barcode_prefix), "", second_barcode_dataframe$Barcode) 
head(second_barcode_dataframe)
nrow(second_barcode_dataframe)

message(sprintf("Creating all barcodes lookup table\n"))
# Get common barcodes between both barcode lists
all_barcodes <- union(first_barcode_dataframe$Barcode, second_barcode_dataframe$Barcode)
#head(all_barcodes)
all_barcodeIDs <- c(1:length(all_barcodes))
#head(all_barcodeIDs)
#http://best-answer.net/easy-way-to-perform-a-lookup-in-r/
all_barcodes_lookup_table <- data.frame("BarcodeID" = all_barcodeIDs, "Barcode" = all_barcodes)
head(all_barcodes_lookup_table)
nrow(all_barcodes_lookup_table)

message(sprintf("Reading first cell by bin summary csv\n"))
first_cell_by_bin_summary_dataframe <- read.csv(opt$first_cell_by_bin_summary_csv)
first_cell_by_bin_summary_dataframe['value'] = 1
# Remove the barcode ID column since we will replace them with barcodes
# that will be common between the two tables
first_cell_by_bin_summary_dataframe$Barcode <- first_barcode_dataframe[match(first_cell_by_bin_summary_dataframe$BarcodeID, first_barcode_dataframe$BarcodeID), "Barcode"]
first_cell_by_bin_summary_dataframe$BarcodeID <- NULL
# Look up the barcode string in the look up table using the barcode ID
# column and add a new column for the barcode string
first_cell_by_bin_summary_dataframe$BarcodeID <- all_barcodes_lookup_table[match(first_cell_by_bin_summary_dataframe$Barcode, all_barcodes_lookup_table$Barcode), "BarcodeID"]
first_cell_by_bin_summary_dataframe$Barcode <- NULL
head(first_cell_by_bin_summary_dataframe)

message(sprintf("Barcode match\n"))
print(match("AACGTGATCTGAGCCAACCTCCAA", all_barcodes_lookup_table$Barcode))
print(match("CCTCTATCAAACATCGACGCTCGA", all_barcodes_lookup_table$Barcode))

message(sprintf("\nBarcodes in first barcode list that are not in the second barcode list\n"))
barcodes_in_first_not_in_second <- as.data.frame(setdiff(first_barcode_dataframe$Barcode, second_barcode_dataframe$Barcode))
barcodes_in_first_not_in_second
nrow(barcodes_in_first_not_in_second)

message(sprintf("\nBarcodes in second barcode list that are not in the first barcode list\n"))
barcodes_in_second_not_in_first <- as.data.frame(setdiff(second_barcode_dataframe$Barcode, first_barcode_dataframe$Barcode))
barcodes_in_second_not_in_first
nrow(barcodes_in_second_not_in_first)

message(sprintf("\nBarcodes in the first barcode list and the second barcode list\n"))
common_barcodes <- as.data.frame(intersect(first_barcode_dataframe$Barcode, second_barcode_dataframe$Barcode))
nrow(common_barcodes)


message(sprintf("\nReading second cell by bin summary csv\n"))
second_cell_by_bin_summary_dataframe <- read.csv(opt$second_cell_by_bin_summary_csv)
second_cell_by_bin_summary_dataframe['value'] = 1
second_cell_by_bin_summary_dataframe$Barcode <- second_barcode_dataframe[match(second_cell_by_bin_summary_dataframe$BarcodeID, second_barcode_dataframe$BarcodeID), "Barcode"]
# Remove the barcode ID column since we will replace them with barcodes
# that will be common between the two tables
second_cell_by_bin_summary_dataframe$BarcodeID <- NULL
# Look up the barcode string in the look up table using the barcode ID
# column and add a new column for the barcode string
second_cell_by_bin_summary_dataframe$BarcodeID <-  all_barcodes_lookup_table[match(second_cell_by_bin_summary_dataframe$Barcode, all_barcodes_lookup_table$Barcode), "BarcodeID"]
head(second_cell_by_bin_summary_dataframe)
second_cell_by_bin_summary_dataframe$Barcode <- NULL
# Look up the barcode string in the look up table using the barcode ID
# column and add a new column for the barcode string
head(second_cell_by_bin_summary_dataframe)

message(sprintf("\nCreating first cell by bin summary matrix\n"))
first_cell_by_bin_summary_matrix <- as.matrix(first_cell_by_bin_summary_dataframe)

head(first_cell_by_bin_summary_matrix)

message(sprintf("\nCreating second cell by bin summary matrix\n"))
second_cell_by_bin_summary_matrix <- as.matrix(second_cell_by_bin_summary_dataframe)

head(second_cell_by_bin_summary_matrix)

# https://stackoverflow.com/questions/47973518/how-to-convert-sparse-matrix-to-dense-matrix-in-r
# create matrix large enough to hold data
#make sure the dense matrix can hold the biggest dimensions of input matrices
max_rows <- max(first_cell_by_bin_summary_dataframe$BarcodeID, second_cell_by_bin_summary_dataframe$BarcodeID)
max_rows
max_columns <- max(first_cell_by_bin_summary_dataframe$Bin, second_cell_by_bin_summary_dataframe$Bin)
max_columns

message(sprintf("\nCreating first cell by bin dense matrix\n"))
first_cell_by_bin_dense_matrix = matrix(0, max_rows, max_columns)

first_cell_by_bin_dense_matrix[cbind(first_cell_by_bin_summary_dataframe$BarcodeID, first_cell_by_bin_summary_dataframe$Bin)] <- first_cell_by_bin_summary_dataframe$value
message(sprintf("\nNumber of first cell by bin dense matrix entries\n"))
sum_first_cell_by_bin_values <- sum(first_cell_by_bin_summary_dataframe$value)
sum_first_cell_by_bin_values
message(sprintf("\nFirst cell by bin dense matrix [30:50,1:10]\n"))
first_cell_by_bin_dense_matrix[30:50,1:10]

message(sprintf("\nCreating second cell by bin dense matrix\n"))
second_cell_by_bin_dense_matrix = matrix(0, max_rows, max_columns)

second_cell_by_bin_dense_matrix[cbind(second_cell_by_bin_summary_dataframe$BarcodeID, second_cell_by_bin_summary_dataframe$Bin)] <- second_cell_by_bin_summary_dataframe$value
message(sprintf("\nNumber of second cell by bin dense matrix entries\n"))
sum_second_cell_by_bin_values <- sum(second_cell_by_bin_summary_dataframe$value)
sum_second_cell_by_bin_values
message(sprintf("\nSecond cell by bin dense matrix [100:120,1:10]\n"))
second_cell_by_bin_dense_matrix[100:120,1:10]


#It would be great to compute the Jaccard index of the sets of barcodes from your version and hers
#(Size of intersection divided by size of union)
#And then the number and proportion of matrix entries that differ in the two cell-by-bin matrices
message(sprintf("\nSum of differences between first and second dense matrices\n"))
sum_of_differences_dense_matrices <- sum(first_cell_by_bin_dense_matrix != second_cell_by_bin_dense_matrix)
sum_of_differences_dense_matrices

message(sprintf("\nProportion of matrix entries that differ between the dense matrices\n"))
proportion_of_differences_dense_matrices <- mean(first_cell_by_bin_dense_matrix != second_cell_by_bin_dense_matrix)
proportion_of_differences_dense_matrices

#message(sprintf("\nSum of intersection between first and second dense matrices\n"))
#interection_first_and_second <- merge(first_cell_by_bin_dense_matrix, second_cell_by_bin_dense_matrix)
#sum_of_intersection_dense_matrices <- sum(interection_first_and_second)
#sum_of_intersection_dense_matrices
#

