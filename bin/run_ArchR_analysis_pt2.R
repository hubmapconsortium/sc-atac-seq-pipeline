#need to reload the libraries
library(optparse)
library(parallel)
library(magick)
library(R.utils)
library(BSgenome.Hsapiens.UCSC.hg38)

option_list <- list(
  make_option(
    c("-i", "--image_file"),
    type = "character",
    help = "Path to the RData image from previous ArchR step"
  ),
  make_option(
    c("-a", "--archr_project"),
    type = "character",
    help = "Path to the ArchR Project directory from previous ArchR step"
  ))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
image_path <- opt$image_file
archr_path <- opt$archr_project
copyDirectory(archr_path, "ArchRStep2")
load(image_path)

library(ArchR)

archr_proj <- loadArchRProject(path = "/output/ArchRStep2")
addArchRGenome("hg38")

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

archr_proj <- saveArchRProject(ArchRProj = archr_proj)

save.image(file="atacSeqStep2.RData")
