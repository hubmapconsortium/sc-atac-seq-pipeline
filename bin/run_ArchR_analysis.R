#!/usr/bin/env Rscript
library(optparse)

# https://github.com/GreenleafLab/ArchR/discussions/
# 1044#discussioncomment-1405648
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

# First, we load the ArchR library. If this fails, you have not properly
# installed ArchR and should revisit the installation instructions. We
# also recommend setting and remembering a known seed to facilitate
# replication of operations requiring randomization.
library(ArchR)

# Next, we set the default number of threads for parallelized operations
# in ArchR functions. You should change the value passed to threads to match
# the specifications of your local machine.
addArchRThreads(threads = opt$threads)

script_dir <- getwd()
message(paste("\n\nDirectory used is:", script_dir, "\n\n"))

input_files <- c(opt$bam_file)
message(paste0("\n\nNames of input files: ", input_files, "\n\n"))

names(input_files) <- c("BAM_data")

# Before we begin, we need add a reference genome annotation for ArchR to have
# access to chromosome and gene information. ArchR natively supports hg19,
# hg38, mm9, and mm10.
addArchRGenome("hg38")

# Creating Arrow Files
# Now we will create our Arrow files. For each sample, this step will:

# Read accessible fragments from the provided input files.
# Calculate quality control information for each cell (i.e. TSS enrichment
# scores and nucleosome info).
# Filter cells based on quality control parameters.
# Create a genome-wide TileMatrix using 500-bp bins.
# Create a GeneScoreMatrix using the custom geneAnnotation that was defined
# when we called addArchRGenome().
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

# Inferring Doublets
# After Arrow file creation, we can infer potential doublets (a single droplet
# containing multiple cells) that can confound downstream results. This is
# done using the addDoubletScores() function.
message(paste0("Adding Doublet Scores"))
archr_proj <- addDoubletScores(
  input = archr_proj,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor
                      # search.
)

# We can check how much memory is used to store the ArchRProject in memory
# within R:
message(paste0("\n\nMemory Size = ",
        round(object.size(archr_proj) / 10^6, 3), " MB\n\n"))
## [1] “Memory Size = 37.135 MB”

# We can also ask which data matrices are available within the ArchRProject
# which will be useful downstream once we start adding to this project:
getAvailableMatrices(archr_proj)
## [1] “GeneScoreMatrix” “TileMatrix”

# We can access the cell names associated with each cell:
head(archr_proj$cellNames)
## [1] “scATAC_BMMC_R1#TTATGTCAGTGATTAG-1” “scATAC_BMMC_R1#AAGATAGTCACCGCGA-1”
## [3] “scATAC_BMMC_R1#GCATTGAAGATTCCGT-1” “scATAC_BMMC_R1#TATGTTCAGGGTTCCC-1”
## [5] “scATAC_BMMC_R1#TCCATCGGTCCCGTGA-1” “scATAC_BMMC_R1#AGTTACGAGAACGTCG-1”

# We can access the sample names associated with each cell:
head(archr_proj$Sample)
## [1] “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1”
## [5] “scATAC_BMMC_R1” “scATAC_BMMC_R1”

#We can access the TSS Enrichment Scores for each cell:
quantile(archr_proj$TSSEnrichment)
## 0% 25% 50% 75% 100%
## 4.027 13.922 16.832 19.937 41.782

# Plot QC metrics - log10(Unique Fragments) vs TSS enrichment score

# Repeating the example shown above, we can easily obtain standard scATAC-seq
# metrics for quality control of individual cells. We have found that the most
# robust metrics for quality control are the TSS enrichment score (a measure of
# signal-to-background in ATAC-seq data) and the number of unique nuclear
# fragments (because cells with very few fragments do not have enough data to
# confidently analyze).
df <- getCellColData(archr_proj, select = c("log10(nFrags)", "TSSEnrichment"))

# Now lets plot the number of unique nuclear fragments (log10) by the TSS
# enrichment score.
# This type of plot is key for identifying high quality cells. You’ll notice
# that the cutoffs that we previously specified when creating the Arrow files
# (via minTSS and minFrags) have already removed low quality cells.
# However, if we noticed that the previously applied QC filters were not
# adequate for this sample, we could further adjust our cutoffs based on this
# plot or re-generate the Arrow files if needed.
p <- ggPoint(
     x = df[, 1],
     y = df[, 2],
     colorDensity = TRUE,
     continuousSet = "sambaNight",
     xlabel = "Log10 Unique Fragments",
     ylabel = "TSS Enrichment",
     xlim = c(log10(500), quantile(df[, 1], probs = 0.99)),
         ylim = c(0, quantile(df[, 2], probs = 0.99))
     ) + geom_hline(yintercept = 4, lty = "dashed") +
           geom_vline(xintercept = 3, lty = "dashed")

#To save an editable vectorized version of this plot, we use plotPDF().
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = archr_proj, addDOC = FALSE)

#Make a ridge plot for each sample for the TSS enrichment scores.
#To make a ridge plot, we set plotAs = "ridges".

p1 <- plotGroups(
     ArchRProj = archr_proj,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "TSSEnrichment",
     plotAs = "ridges"
        )

# Make a violin plot for each sample for the TSS enrichment scores.
# To make a violin plot, we set plotAs = "violin". Violin plots in ArchR come
# with a box-and-whiskers plot in the style of Tukey as implemented by ggplot2.
# This means that the lower and upper hinges correspond to the 25th and 75th
# percentiles, respectively, and the middle corresponds to the median.i
# The lower and upper whiskers extend from the hinge to the lowest or highest
# value or 1.5 times the interquartile range (the distance between the 25th
# and 75th percentiles).
p2 <- plotGroups(
     ArchRProj = archr_proj,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "TSSEnrichment",
     plotAs = "violin",
         alpha = 0.4,
         addBoxPlot = TRUE
    )

# Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
     ArchRProj = archr_proj,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "log10(nFrags)",
     plotAs = "ridges"
        )

# Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
     ArchRProj = archr_proj,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "log10(nFrags)",
     plotAs = "violin",
     alpha = 0.4,
     addBoxPlot = TRUE
    )

# To save editable vectorized versions of these plots, we use plotPDF().
plotPDF(p1, p2, p3, p4, name = "QC-Sample-Statistics.pdf",
       ArchRProj = archr_proj, addDOC = FALSE, width = 4, height = 4)

# Plot Sample Fragment Size Distribution and TSS Enrichment Profiles.
# Because of how the data is stored and accessed, ArchR can compute fragment
# size distributions and TSS enrichment profiles from Arrow files very quickly.

# To plot the fragment size distributions of all samples, we use the
# plotFragmentSizes() function. Fragment size distributions in ATAC-seq can be
# quite variable across samples, cell types, and batches. Slight differences
# like those shown below are common and do not necessarily correlate with
# differences in data quality.
pfrag <- plotFragmentSizes(ArchRProj = archr_proj)

# Plot TSS enrichment profiles, We use the plotTSSEnrichment() function. TSS
# enrichment
# profiles should show a clear peak in the center and a smaller shoulder peak
# right-of-center which is caused by the well-positioned +1 nucleosome.
ptssen <- plotTSSEnrichment(ArchRProj = archr_proj)

#To save editable vectorized versions of these plots, we use plotPDF().
plotPDF(pfrag, ptssen, name = "QC-Sample-FragSizes-TSSProfile.pdf",
       ArchRProj = archr_proj, addDOC = FALSE, width = 5, height = 5)

# Now we can filter putative doublets based on the previously determined
# doublet scores using the filterDoublets() function. This doesn’t physically
# remove data from the Arrow files but rather tells the ArchRProject to ignore
# these cells for downstream analysis.
message(paste("Filtering Doublets"))
archr_proj <- filterDoublets(ArchRProj = archr_proj)

## Dimensionality Reduction and Clustering
## ArchR implements an iterative LSI dimensionality reduction via the
# addIterativeLSI() function.
message(paste("Doing Dimensionality Reduction"))
archr_proj <- addIterativeLSI(
    ArchRProj = archr_proj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 5,
    clusterParams = list(
        resolution = c(2),
        sampleCells = 10000,
        maxClusters = 6,
        n.start = 10),
    varFeatures = 25000
    )

# To call clusters in this reduced dimension sub-space, we use the addClusters()
# function which uses Seurat’s graph clustering as the default clustering
# method.
message(paste("Adding Clusters"))
archr_proj <- addClusters(input = archr_proj, reducedDims = "IterativeLSI")

# Visualizing in a 2D UMAP Embedding
# We can visualize our scATAC-seq data using a 2-dimensional representation
# such as Uniform Manifold Approximation and Projection (UMAP). To do this, we
# add a UMAP embedding to our ArchRProject object with the addUMAP() function.
# This function uses the uwot package to perform UMAP.

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


message(paste("Adding UMAP"))
archr_proj <- addUMAP(ArchRProj = archr_proj, reducedDims = "IterativeLSI")

message(paste("Getting embedding"))
archr_proj_embed_w_clusters_df <- getEmbedding(ArchRProj = archr_proj,
       embedding = "UMAP", returnDF = TRUE)

message(paste("Adding Clusters column"))
# https://stackoverflow.com/questions/48896190/
# add-column-to-r-dataframe-based-on-rowname
# https://intellipaat.com/community/31833/
# r-add-a-new-column-to-a-dataframe-using-matching-values-of-another-dataframe
# row.names(archr_proj_embed_w_clusters_df) must be the first argument because
# cell_col_data_df may have rows that do not exist in
# archr_proj_embed_w_clusters_df and match in that case will return a larger
# vector with NAs, which will fail when used as an index to
# cell_col_data_df$Clusters.
archr_proj_embed_w_clusters_df$Clusters <- cell_col_data_df$Clusters[
     match(row.names(archr_proj_embed_w_clusters_df),
    row.names(cell_col_data_df))]
write.csv(archr_proj_embed_w_clusters_df,
          file = "umap_coords_clusters.csv")


# Using this UMAP, we can visualize various attributes of our cells which are
# stored in a matrix called cellColData in our ArchRProject. To do this, we use
# the plotEmbedding() function and we specify the variable to use for coloration
# via a combination of the colorBy and name parameters.
#
## For example, we can color by “Sample”:
pcell_col_data_sample_umap <- plotEmbedding(ArchRProj = archr_proj,
            colorBy = "cellColData", name = "Sample", embedding = "UMAP")

# Or we can color by “Clusters”:
pcell_col_data_clusters_umap <- plotEmbedding(ArchRProj = archr_proj,
             colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(pcell_col_data_sample_umap, pcell_col_data_clusters_umap,
             type = "h")
# To save an editable vectorized version of this plot, we use the plotPDF()
# function.
plotPDF(pcell_col_data_sample_umap, pcell_col_data_clusters_umap,
        name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = archr_proj, addDOC = FALSE, width = 5, height = 5)

## Assigning Clusters with Gene Scores
# First, we add imputation weights using MAGIC to help smooth the dropout noise
# in our gene scores.
archr_proj <- addImputeWeights(archr_proj)

# Create hdf5 file containing smoothed data
message(paste("Writing smoothed data to hdf5 file"))
smoothed_gene_score_matrix_se <- getMatrixFromProject(ArchRProj = archr_proj,
                                     useMatrix = "GeneScoreMatrix")
smoothed_gene_score_dm <- assays(smoothed_gene_score_matrix_se)$GeneScoreMatrix
smooth_cell_by_gene_filename <- "cell_by_gene_smoothed.hdf5"
message(paste("Writing smoothed cell by gene data to",
                smooth_cell_by_gene_filename))
h5createFile(smooth_cell_by_gene_filename)
# AnnData expects barcodes as rows not columns in convert_to_h5ad.cwl
transposed_smooth_g_score_dm <- t(smoothed_gene_score_dm)
h5write(as.matrix(transposed_smooth_g_score_dm), smooth_cell_by_gene_filename,
                  "cell_by_gene_smoothed", level = 0)
h5write(gene_row_dat_df$name, smooth_cell_by_gene_filename, "genes")
h5write(archr_proj$cellNames, smooth_cell_by_gene_filename, "barcodes")


archr_proj <- addGroupCoverages(ArchRProj = archr_proj, groupBy = "Clusters")
path_to_macs2 <- findMacs2()
archr_proj <- addReproduciblePeakSet(
    ArchRProj = archr_proj,
    groupBy = "Clusters",
    pathToMacs2 = path_to_macs2
    )

peaks_gr <- getPeakSet(archr_proj)
message(paste("Writing peaks CSV and BED files"))
write.csv(peaks_gr, file = "peaks.csv")
library(rtracklayer)
export.bed(peaks_gr, con = "peaks.bed")

archr_proj <- addPeakMatrix(archr_proj)

message(paste("Cell types:"))
# First, lets remind ourselves of the cell types that we are working with in the
# project and their relative proportions.
table(archr_proj$Clusters)

archr_proj <- saveArchRProject(ArchRProj = archr_proj)

# To identify marker genes based on gene scores, we call the
# getMarkerFeatures() function with useMatrix = "GeneScoreMatrix".
# We specify that we want to know the cluster-specific features
# with groupBy = "Clusters" which tells ArchR to use
# the “Clusters” column in cellColData to stratify cell groups.
marker_gs <- getMarkerFeatures(
       ArchRProj = archr_proj,
       useMatrix = "GeneScoreMatrix",
       groupBy = "Clusters",
       bias = c("TSSEnrichment", "log10(nFrags)"),
       testMethod = "wilcoxon"
 )

markers_gs_list <- getMarkers(marker_gs,
             cutOff = "FDR <= 0.01 & Log2FC >= .5")
markers_gs_list

message(paste("Writing gene markers CSV"))
write.csv(markers_gs_list, file = "gene_markers.csv")
  
if (!isEmpty(markers_gs_list)) {
    heatmap_gs <- plotMarkerHeatmap(
         seMarker = marker_gs,
         cutOff = "FDR <= 0.01 & Log2FC >= .5",
         transpose = TRUE
      )
    draw(heatmap_gs, heatmap_legend_side = "bot",
            annotation_legend_side = "bot")
    plotPDF(heatmap_gs, name = "GeneScores-Marker-Heatmap", width = 8,
              height = 6, ArchRProj = archr_proj, addDOC = FALSE)
} else {
    message("No markers found so no marker gene scores heatmap can be created.")
}

# Often times, we are interested to know which peaks are unique to an individual
# cluster or a small group of clusters. We can do this in an unsupervised
# fashion in ArchR using the addMarkerFeatures() function in combination with
# useMatrix = "PeakMatrix". Now, we are ready to identify marker peaks by
# calling the addMarkerFeatures() function with useMatrix = "PeakMatrix".
# Additionally, we tell ArchR to account for differences in data quality amongst
# the cell groups by setting the bias parameter to account for TSS enrichment
# and the number of unique fragments per cell.
markers_peaks <- getMarkerFeatures(
    ArchRProj = archr_proj,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
    )

# The object returned by the getMarkerFeatures() function is a
# SummarizedExperiment that contains a few different assays.
markers_peaks

# Instead of a list of DataFrame objects, we can use getMarkers() to return a
# GRangesList object by setting returnGR = TRUE.
markers_gr <- getMarkers(markers_peaks, cutOff = "FDR <= 0.01 & Log2FC >= .5",
                          returnGR = TRUE)
markers_gr
write.csv(markers_gr, file = "peak_markers.csv")

# ArchR provides multiple plotting functions to interact with the
# SummarizedExperiment objects returned by getMarkerFeatures().
# We can visualize these marker peaks (or any features output by
# getMarkerFeatures()) as a heatmap using the markerHeatmap() function.
if (!isEmpty(markers_gr)) {
    heatmap_peaks <- plotMarkerHeatmap(
      seMarker = markers_peaks,
      cutOff = "FDR <= 0.01 & Log2FC >= .5",
      transpose = TRUE
      )
    # We can plot this heatmap using draw().
    draw(heatmap_peaks, heatmap_legend_side = "bot",
                              annotation_legend_side = "bot")
    plotPDF(heatmap_peaks, name = "Peak-Marker-Heatmap", width = 8, height = 6,
            ArchRProj = archr_proj, addDOC = FALSE)
} else {
    message("No markers found so no marker peaks heatmap can be created.")
}

# Saving and Loading an ArchRProject
# To easily save an ArchRProject for later use or for sharing with
# collaborators,we use the saveArchRProject() function. This copies the current
# ArchRProject object and all of the Arrow files to a specified directory. If
# we don’t specify an output directory (as below), saveArchRProject() uses
# the output directory that we specified upon creation of our ArchRProject.
# In this case that is the folder "ArchRProjFiles"
archr_proj <- saveArchRProject(ArchRProj = archr_proj)

# Session Information
# This tutorial was run on the date specified below.

Sys.Date()

# The sessionInfo() at run time was:
sessionInfo()
