#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(
    c("-b", "--bam_file"),
    type="character",
    default=NULL,
    help="Selected barcodes"
  ),
  make_option(
    c("-n", "--threads"),
    type="integer",
    default=2,
    help="Number of subprocesses/threads to use"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$bam_file)){
  print_help(opt_parser)
  stop("--bam_file argument must be supplied (input BAM file).", call.=FALSE)
}

# First, we load the ArchR library. If this fails, you have not properly installed
# ArchR and should revisit the installation instructions. We also recommend setting
# and remembering a known seed to facilitate replication of operations requiring randomization.
library(ArchR)

# Next, we set the default number of threads for parallelized operations in ArchR
# functions. You should change the value passed to threads to match the
# specifications of your local machine.
addArchRThreads(threads = opt$threads) 

script.dir <- getwd()
message(paste("directory used is:",script.dir))

inputFiles <- c(opt$bam_file)
names(inputFiles) <- c("BAM_data")
inputFiles

# Before we begin, we need add a reference genome annotation for ArchR to have
# access to chromosome and gene information. ArchR natively supports hg19, hg38, mm9, and mm10.
addArchRGenome("hg38")

# Creating Arrow Files
# Now we will create our Arrow files. For each sample, this step will:

# Read accessible fragments from the provided input files.
# Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
# Filter cells based on quality control parameters.
# Create a genome-wide TileMatrix using 500-bp bins.
# Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 1.5, #Dont set this too high because you can always increase later
  minFrags = 2000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  bamFlag = list(isMinusStrand = FALSE, isProperPair = TRUE, isDuplicate = FALSE),
  bcTag = "CB" # We added this tag to the SAM file and then converted it to a BAM
)

projSci <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "ArchRProjFiles",
  copyArrows = TRUE # This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# We can check how much memory is used to store the ArchRProject in memory within R:
paste0("Memory Size = ", round(object.size(projSci) / 10^6, 3), " MB")
## [1] “Memory Size = 37.135 MB”

# We can also ask which data matrices are available within the ArchRProject which
# will be useful downstream once we start adding to this project:
getAvailableMatrices(projSci)
## [1] “GeneScoreMatrix” “TileMatrix”


# We can access the cell names associated with each cell:
head(projSci$cellNames)
## [1] “scATAC_BMMC_R1#TTATGTCAGTGATTAG-1” “scATAC_BMMC_R1#AAGATAGTCACCGCGA-1”
## [3] “scATAC_BMMC_R1#GCATTGAAGATTCCGT-1” “scATAC_BMMC_R1#TATGTTCAGGGTTCCC-1”
## [5] “scATAC_BMMC_R1#TCCATCGGTCCCGTGA-1” “scATAC_BMMC_R1#AGTTACGAGAACGTCG-1”

# We can access the sample names associated with each cell:
head(projSci$Sample)
## [1] “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1”
## [5] “scATAC_BMMC_R1” “scATAC_BMMC_R1”

#We can access the TSS Enrichment Scores for each cell:
quantile(projSci$TSSEnrichment)
## 0% 25% 50% 75% 100%
## 4.027 13.922 16.832 19.937 41.782

# Plot QC metrics - log10(Unique Fragments) vs TSS enrichment score

# Repeating the example shown above, we can easily obtain standard scATAC-seq
# metrics for quality control of individual cells. We have found that the most
# robust metrics for quality control are the TSS enrichment score (a measure of
# signal-to-background in ATAC-seq data) and the number of unique nuclear fragments
# (because cells with very few fragments do not have enough data to confidently analyze).
df <- getCellColData(projSci, select = c("log10(nFrags)", "TSSEnrichment"))

# Now lets plot the number of unique nuclear fragments (log10) by the TSS enrichment score. 
# This type of plot is key for identifying high quality cells. You’ll notice
# that the cutoffs that we previously specified when creating the Arrow files
# (via minTSS and minFrags) have already removed low quality cells. 
# However, if we noticed that the previously applied QC filters were not adequate
# for this sample, we could further adjust our cutoffs based on this plot or
# re-generate the Arrow files if needed.
p <- ggPoint(
     x = df[,1],
     y = df[,2],
     colorDensity = TRUE,
     continuousSet = "sambaNight",
     xlabel = "Log10 Unique Fragments",
     ylabel = "TSS Enrichment",
     xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
         ylim = c(0, quantile(df[,2], probs = 0.99))
     ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

#To save an editable vectorized version of this plot, we use plotPDF().
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projSci, addDOC = FALSE)

#Make a ridge plot for each sample for the TSS enrichment scores.
#To make a ridge plot, we set plotAs = "ridges".

p1 <- plotGroups(
     ArchRProj = projSci,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "TSSEnrichment",
     plotAs = "ridges"
        )

# Make a violin plot for each sample for the TSS enrichment scores.
# To make a violin plot, we set plotAs = "violin". Violin plots in ArchR come
# with a box-and-whiskers plot in the style of Tukey as implemented by ggplot2.
# This means that the lower and upper hinges correspond to the 25th and 75th
# percentiles, respectively, and the middle corresponds to the median. The lower
# and upper whiskers extend from the hinge to the lowest or highest value or 1.5
# times the interquartile range (the distance between the 25th and 75th percentiles).
p2 <- plotGroups(
     ArchRProj = projSci, 
     groupBy = "Sample", 
     colorBy = "cellColData", 
     name = "TSSEnrichment",
     plotAs = "violin",
         alpha = 0.4,
         addBoxPlot = TRUE
    )

# Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
     ArchRProj = projSci, 
     groupBy = "Sample", 
     colorBy = "cellColData", 
     name = "log10(nFrags)",
     plotAs = "ridges"
        )

# Make a violin plot for each sample for the log10(unique nuclear fragments).
p4 <- plotGroups(
     ArchRProj = projSci, 
     groupBy = "Sample", 
     colorBy = "cellColData", 
     name = "log10(nFrags)",
     plotAs = "violin",
     alpha = 0.4,
     addBoxPlot = TRUE
    )

# To save editable vectorized versions of these plots, we use plotPDF().
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projSci, addDOC = FALSE, width = 4, height = 4)

# Plot Sample Fragment Size Distribution and TSS Enrichment Profiles.
# Because of how the data is stored and accessed, ArchR can compute fragment
# size distributions and TSS enrichment profiles from Arrow files very quickly.

# To plot the fragment size distributions of all samples, we use the plotFragmentSizes()
# function. Fragment size distributions in ATAC-seq can be quite variable across
# samples, cell types, and batches. Slight differences like those shown below
# are common and do not necessarily correlate with differences in data quality.
pfrag <- plotFragmentSizes(ArchRProj = projSci)

# Plot TSS enrichment profiles, We use the plotTSSEnrichment() function. TSS enrichment
# profiles should show a clear peak in the center and a smaller shoulder peak
# right-of-center which is caused by the well-positioned +1 nucleosome.
pTSSEn <- plotTSSEnrichment(ArchRProj = projSci)

#To save editable vectorized versions of these plots, we use plotPDF().
plotPDF(pfrag,pTSSEn, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projSci, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projSci, outputDirectory = "ArchRProjFiles", load = FALSE)

# Inferring Doublets
# After Arrow file creation, we can infer potential doublets (a single droplet
# containing multiple cells) that can confound downstream results. This is
# done using the addDoubletScores() function.
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  dimsToUse = 1:15 # Have to make upper dimension less than default of 30 since we only have 18 columns... 
)

# We can also ask which data matrices are available within the ArchRProject
# which will be useful downstream once we start adding to this project:
getAvailableMatrices(projSci)
### [1] “GeneScoreMatrix” “TileMatrix”

message(paste("Tile Matrix:"))
tile_matrix = getMatrixFromProject(
     ArchRProj = projSci,
     useMatrix = "TileMatrix",
     useSeqnames = "chr1",
     verbose = TRUE,
     binarize = TRUE,
     threads = getArchRThreads(),
     logFile = createLogFile("getTileMatrixFromProject")
   )
tile_matrix

message(paste("Gene Score Matrix:"))
gene_score_matrix = getMatrixFromProject(
     ArchRProj = projSci,
     useMatrix = "GeneScoreMatrix",
     useSeqnames = "chr1",
     verbose = TRUE,
     binarize = TRUE,
     threads = getArchRThreads(),
     logFile = createLogFile("getGeneScoreMatrixFromProject")
   )
gene_score_matrix


## Now we can filter putative doublets based on the previously determined
## doublet scores using the filterDoublets() function. This doesn’t physically
## remove data from the Arrow files but rather tells the ArchRProject to ignore
# these cells for downstream analysis.
#
projSci <- filterDoublets(ArchRProj = projSci)

## Dimensionality Reduction and Clustering
## ArchR implements an iterative LSI dimensionality reduction via the addIterativeLSI() function.
projSci <- addIterativeLSI(ArchRProj = projSci, useMatrix = "TileMatrix", name = "IterativeLSI")
#
# To call clusters in this reduced dimension sub-space, we use the addClusters()
# function which uses Seurat’s graph clustering as the default clustering method.
projSci <- addClusters(input = projSci, reducedDims = "IterativeLSI")

# Visualizing in a 2D UMAP Embedding
# We can visualize our scATAC-seq data using a 2-dimensional representation
# such as Uniform Manifold Approximation and Projection (UMAP). To do this, we
# add a UMAP embedding to our ArchRProject object with the addUMAP() function.
# This function uses the uwot package to perform UMAP.

cellColDataDF <- getCellColData(projSci)
write.csv(cellColDataDF, file='cell_column_data.csv')

projSci <- addUMAP(ArchRProj = projSci, reducedDims = "IterativeLSI")
projSciEmbeddingWClustersDF = getEmbedding(ArchRProj = projSci, embedding = "UMAP", returnDF = TRUE)
projSciEmbeddingWClustersDF$Clusters <- cellColDataDF$Clusters[match(row.names(cellColDataDF), row.names(projSciEmbeddingWClustersDF))]
write.csv(projSciEmbeddingWClustersDF, file='archr_umap_coords_clusters.csv')


# Using this UMAP, we can visualize various attributes of our cells which are
# stored in a matrix called cellColData in our ArchRProject. To do this, we use
# the plotEmbedding() function and we specify the variable to use for coloration
#via a combination of the colorBy and name parameters.
#
## For example, we can color by “Sample”:
pcellCollDataSampleUMAP <- plotEmbedding(ArchRProj = projSci, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

# Or we can color by “Clusters”:
pCellCollDataClustersUMAP <- plotEmbedding(ArchRProj = projSci, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

ggAlignPlots(pcellCollDataSampleUMAP, pCellCollDataClustersUMAP, type = "h")
# To save an editable vectorized version of this plot, we use the plotPDF() function.
plotPDF(pcellCollDataSampleUMAP, pCellCollDataClustersUMAP, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = projSci, addDOC = FALSE, width = 5, height = 5)

## Assigning Clusters with Gene Scores
# First, we add imputation weights using MAGIC to help smooth the dropout noise in our gene scores.
projSci <- addImputeWeights(projSci)
### ArchR logging to : ArchRLogs/ArchR-addImputeWeights-69ef433c71d0-Date-2020-04-21_Time-16-33-19.log
### If there is an issue, please report to github with logFile!
### 2020-04-21 16:33:19 : Computing Impute Weights Using Magic (Cell 2018), 0 mins elapsed.

projSci <- addGroupCoverages(ArchRProj = projSci, groupBy = "Clusters")
pathToMacs2 <- findMacs2()
projSci <- addReproduciblePeakSet(
    ArchRProj = projSci, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
    )

peaks_gr = getPeakSet(projSci)
write.csv(peaks_gr, file = "peaks.csv")

projSci <- addPeakMatrix(projSci)
getAvailableMatrices(projSci)


message(paste("Cell types:"))
# First, lets remind ourselves of the cell types that we are working with in the
# project and their relative proportions.
table(projSci$Clusters)


# To identify marker genes based on gene scores, we call the getMarkerFeatures()
# function with useMatrix = "GeneScoreMatrix". We specify that we want to know
# the cluster-specific features with groupBy = "Clusters" which tells ArchR to use
# the “Clusters” column in cellColData to stratify cell groups.
markersGS <- getMarkerFeatures(
       ArchRProj = projSci, 
       useMatrix = "GeneScoreMatrix", 
       groupBy = "Clusters",
       bias = c("TSSEnrichment", "log10(nFrags)"),
       testMethod = "wilcoxon"
 )

markerGSList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
#markerGSList$C1

heatmapGS <- markerHeatmap(
     seMarker = markersGS, 
     cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
     transpose = TRUE
  )

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projSci, addDOC = FALSE)


# Often times, we are interested to know which peaks are unique to an individual
# cluster or a small group of clusters. 
# We can do this in an unsupervised fashion in ArchR using the addMarkerFeatures()
# function in combination with useMatrix = "PeakMatrix".
# Now, we are ready to identify marker peaks by calling the addMarkerFeatures() function with useMatrix = "PeakMatrix". 
# Additionally, we tell ArchR to account for differences in data quality amongst the cell groups by setting the bias
# parameter to account for TSS enrichment and the number of unique fragments per cell.
markersPeaks <- getMarkerFeatures(
    ArchRProj = projSci, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
    )

# The object returned by the getMarkerFeatures() function is a SummarizedExperiment that contains a few different assays.
markersPeaks

# Instead of a list of DataFrame objects, we can use getMarkers() to return a
# GRangesList object by setting returnGR = TRUE.
markers_gr <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markers_gr
write.csv(markers_gr, file = "markers.csv")

# ArchR provides multiple plotting functions to interact with the SummarizedExperiment objects returned by getMarkerFeatures().
# We can visualize these marker peaks (or any features output by getMarkerFeatures()) as a heatmap using the markerHeatmap() function.
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
  )


# We can plot this heatmap using draw().
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projSci, addDOC = FALSE)

# Saving and Loading an ArchRProject
# To easily save an ArchRProject for later use or for sharing with collaborators,
# we use the saveArchRProject() function. This copies the current ArchRProject object
# and all of the Arrow files to a specified directory. If we don’t specify an output
# directory (as below), saveArchRProject() uses the output directory that we specified
# upon creation of our ArchRProject. In this case that is the folder 'ArchRProjFiles'
projSci <- saveArchRProject(ArchRProj = projSci)

# When we are ready to load this saved ArchRProject we use the loadArchRProject()
# object and provide the path to the folder containing the saved ArchRProject object.
#projSci <- loadArchRProject(path = "ArchRProjFiles")

# Session Information
# This tutorial was run on the date specified below.

Sys.Date()
## [1] “2020-04-21”

# The sessionInfo() at run time was:
sessionInfo()

