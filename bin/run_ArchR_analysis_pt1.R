#!/usr/bin/env Rscript
library(optparse)

# Following two library calls needed when using R 4.1.1
library(parallel)
library(magick)

option_list <- list(
  make_option(
    c("-b", "--bam_file"),
    type = "character",
    default = NULL,
    help = "BAM file to use."
  ),
  make_option(
    c("-i", "--bam_index"),
    type="character",
    default=NULL,
    help="BAM file index."
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 2,
    help = "Number of subprocesses/threads to use."
  ),
  make_option(
    c("-e", "--minTSS"),
    type = "double",
    default = 1.5,
    help <- paste("The minimum numeric transcription start site (TSS)",
                  " enrichment score required to pass filtering. E.g. 1.5", sep = "")
  ),
  make_option(
    c("-g", "--minFrags"),
    type = "integer",
    default = 2000,
    help <- paste("The minimum number of mapped ATAC-seq fragments required",
                  "per cell to pass filtering. E.g. 2000", sep = "")
  ),
  make_option(
    c("-c", "--minCells"),
    type = "integer",
    default = 1000,
    help <- paste("The minimum number of cells in the ArchR project that must",
                  " pass filtering before a warning message is printed. E.g. 1000", sep = "")
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bam_file)) {
  print_help(opt_parser)
  stop("--bam_file argument must be supplied (input BAM file).", call. = FALSE)
}

library(ArchR)

addArchRThreads(threads = opt$threads)

script_dir <- getwd()
message(paste("\n\nDirectory used is:", script_dir, "\n\n"))

input_files <- c(opt$bam_file)
message(paste0("\n\nNames of input files: ", input_files, "\n\n"))

names(input_files) <- c("BAM_data")

addArchRGenome("hg38")

# Create Arrow Files
arrow_files <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  minTSS = opt$minTSS, # Dont set this too high because you can always
  # increase later
  minFrags = opt$minFrags,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  bamFlag = list(isMinusStrand = FALSE, isProperPair = TRUE,
                 isDuplicate = FALSE),
  bcTag = "CB" # We added this tag to the SAM file and then converted
  # it to a BAM
)
arrow_files

archr_proj <- ArchRProject(
  ArrowFiles = arrow_files,
  outputDirectory = "ArchRProjFiles",
  copyArrows = TRUE # This is recommended so that if you modify the Arrow files
  # you have an original copy for later usage.
)
archr_proj

num_cells_pass_filter <- nCells(archr_proj)
message(paste0("\nNumber of cells in the project that passed filtering = ",
               num_cells_pass_filter, "\n\n"))
if (num_cells_pass_filter < opt$minCells) {
  message(paste0("\nWARNING: THE NUMBER OF CELLS IN THE PROJECT IS",
                 " LESS THAN ", opt$minCells,
                 "; THE PIPELINE MAY FAIL UNEXPECTEDLY!\n\n"))
}

# Infer Doublets
message(paste0("Adding Doublet Scores"))
archr_proj <- addDoubletScores(
  input = archr_proj,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor
  # search.
)

message(paste0("\n\nMemory Size = ",
               round(object.size(archr_proj) / 10^6, 3), " MB\n\n"))

getAvailableMatrices(archr_proj)

head(archr_proj$cellNames)

head(archr_proj$Sample)

quantile(archr_proj$TSSEnrichment)

# Plot QC metrics - log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(archr_proj, select = c("log10(nFrags)", "TSSEnrichment"))

#either need to actually plot these or take out of the cwl

# Filter doublets
message(paste("Filtering Doublets"))
archr_proj <- filterDoublets(ArchRProj = archr_proj)

cell_col_data_df <- getCellColData(archr_proj)
write.csv(cell_col_data_df, file = "cell_column_data.csv")

## Create the cell by gene table MTX and CSVs
message(paste("Creating cell by gene MTX file"))
gene_score_matrix_se <- getMatrixFromProject(
  ArchRProj = archr_proj,
  useMatrix = "GeneScoreMatrix",
  logFile = createLogFile("getGeneScoreMatrixFromProject")
)
gene_score_dm <- assays(gene_score_matrix_se)$GeneScoreMatrix
# AnnData expects barcodes as rows not columns in convert_to_h5ad.cwl
transposed_gene_score_dm <- t(gene_score_dm)
writeMM(transposed_gene_score_dm, "cell_by_gene_raw.mtx")

message(paste("Creating gene row data CSV file"))
gene_row_dat_df <- rowData(gene_score_matrix_se)
write.csv(gene_row_dat_df, file = "gene_row_data.csv")

## Create the cell by bin table MTX and CSVs
message(paste("Creating cell by bin MTX file"))
tile_matrix_se <- getMatrixFromProject(
  ArchRProj = archr_proj,
  useMatrix = "TileMatrix",
  binarize = TRUE,
  logFile = createLogFile("getTileMatrixFromProject")
)
tile_dm <- assays(tile_matrix_se)$TileMatrix
# AnnData expects barcodes as rows not columns in convert_to_h5ad.cwl
transposed_tile_dm <- t(tile_dm)
writeMM(transposed_tile_dm, "cell_by_bin.mtx")

message(paste("Creating cell by bin column data CSV file"))
tile_col_data_df <- colData(tile_matrix_se)
write.csv(tile_col_data_df, file = "cell_by_bin_col_data.csv")

message(paste("Creating cell by bin row data CSV file"))
tile_row_data_df <- rowData(tile_matrix_se)
write.csv(tile_row_data_df, file = "cell_by_bin_row_data.csv")

write.table(archr_proj$cellNames, "barcodes.txt", col.names = FALSE,
            row.names = FALSE, quote = FALSE)
write.table(tile_row_data_df, "bins.txt", col.names = FALSE, row.names = FALSE,
            quote = FALSE)


# Create hdf5 file containing smoothed data
message(paste("Writing smoothed data to hdf5 file"))
smoothed_gene_score_matrix_se <- getMatrixFromProject(ArchRProj = archr_proj,
                                                      useMatrix = "GeneScoreMatrix")
smoothed_gene_score_dm <- assays(smoothed_gene_score_matrix_se)$GeneScoreMatrix
smooth_cell_by_gene_filename <- "cell_by_gene_smoothed.hdf5"
message(paste("Writing smoothed cell by gene data to",
              smooth_cell_by_gene_filename))
if(!file.exists(smooth_cell_by_gene_filename)){
  h5createFile(smooth_cell_by_gene_filename)
}
# AnnData expects barcodes as rows not columns in convert_to_h5ad.cwl
transposed_smooth_g_score_dm <- t(smoothed_gene_score_dm)
h5write(as.matrix(transposed_smooth_g_score_dm), smooth_cell_by_gene_filename,
        "cell_by_gene_smoothed", level = 0)
h5write(gene_row_dat_df$name, smooth_cell_by_gene_filename, "genes")
h5write(archr_proj$cellNames, smooth_cell_by_gene_filename, "barcodes")

save.image(file="atacSeqStep1.RData")
saveArchRProject(ArchRProj = archr_proj)
