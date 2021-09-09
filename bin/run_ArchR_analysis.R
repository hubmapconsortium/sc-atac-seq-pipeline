#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(
    c("-b", "--bam_file"),
    type="character",
    default=NULL,
    help="Selected barcodes"
  ),
#  make_option(
#    c("-q", "--QCDir"),
#    type="character",
#    default="QualityControl",
#    help="The relative path to the output directory for QC-level information and plots for each sample/ArrowFile."
#  ),
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
set.seed(1)

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
# Now we will create our Arrow files which will take 10-15 minutes. For each sample, this step will:

# Read accessible fragments from the provided input files.
# Calculate quality control information for each cell (i.e. TSS enrichment scores and nucleosome info).
# Filter cells based on quality control parameters.
# Create a genome-wide TileMatrix using 500-bp bins.
# Create a GeneScoreMatrix using the custom geneAnnotation that was defined when we called addArchRGenome().
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 2000,
#  QCDir = opt$QCDir,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  bamFlag = list(isMinusStrand = FALSE, isProperPair = TRUE, isDuplicate = FALSE),
  bcTag = "CB" # We added this tag to the SAM file and then converted it to a BAM
)

# We can inspect the ArrowFiles object to see that it is actually just a character vector of Arrow file paths.
ArrowFiles
## “scATAC_BMMC_R1.arrow” “scATAC_CD34_BMMC_R1.arrow”
## “scATAC_PBMC_R1.arrow”

projSci <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "SciTest",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

#We can examine the contents of our ArchRProject:
projSci

#We can check how much memory is used to store the ArchRProject in memory within R:
paste0("Memory Size = ", round(object.size(projSci) / 10^6, 3), " MB")
## [1] “Memory Size = 37.135 MB”

#We can also ask which data matrices are available within the ArchRProject which will be useful downstream once we start adding to this project:
getAvailableMatrices(projSci)
## [1] “GeneScoreMatrix” “TileMatrix”


#We can access the cell names associated with each cell:
head(projSci$cellNames)
## [1] “scATAC_BMMC_R1#TTATGTCAGTGATTAG-1” “scATAC_BMMC_R1#AAGATAGTCACCGCGA-1”
## [3] “scATAC_BMMC_R1#GCATTGAAGATTCCGT-1” “scATAC_BMMC_R1#TATGTTCAGGGTTCCC-1”
## [5] “scATAC_BMMC_R1#TCCATCGGTCCCGTGA-1” “scATAC_BMMC_R1#AGTTACGAGAACGTCG-1”

#We can access the sample names associated with each cell:
head(projSci$Sample)
## [1] “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1” “scATAC_BMMC_R1”
## [5] “scATAC_BMMC_R1” “scATAC_BMMC_R1”

#We can access the TSS Enrichment Scores for each cell:
quantile(projSci$TSSEnrichment)
## 0% 25% 50% 75% 100%
## 4.027 13.922 16.832 19.937 41.782

#Example 5. Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score

#Repeating the example shown above, we can easily obtain standard scATAC-seq
# metrics for quality control of individual cells. We have found that the most
# robust metrics for quality control are the TSS enrichment score (a measure of
# signal-to-background in ATAC-seq data) and the number of unique nuclear fragments
# (because cells with very few fragments do not have enough data to confidently analyze).
df <- getCellColData(projSci, select = c("log10(nFrags)", "TSSEnrichment"))
df

#Now lets plot the number of unique nuclear fragments (log10) by the TSS enrichment score. This type of plot is key for identifying high quality cells. You’ll notice that the cutoffs that we previously specified when creating the Arrow files (via filterTSS and filterFrags) have already removed low quality cells. However, if we noticed that the previously applied QC filters were not adequate for this sample, we could further adjust our cutoffs based on this plot or re-generate the Arrow files if needed.
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

p

#To save an editable vectorized version of this plot, we use plotPDF().
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projSci, addDOC = FALSE)

#Example 1. Make a ridge plot for each sample for the TSS enrichment scores.
#To make a ridge plot, we set plotAs = "ridges".

p1 <- plotGroups(
     ArchRProj = projSci,
     groupBy = "Sample",
     colorBy = "cellColData",
     name = "TSSEnrichment",
     plotAs = "ridges"
        )
p1

# Make a violin plot for each sample for the TSS enrichment scores.
# To make a violin plot, we set plotAs = "violin". Violin plots in ArchR come with a box-and-whiskers plot in the style of Tukey as implemented by ggplot2. This means that the lower and upper hinges correspond to the 25th and 75th percentiles, respectively, and the middle corresponds to the median. The lower and upper whiskers extend from the hinge to the lowest or highest value or 1.5 times the interquartile range (the distance between the 25th and 75th percentiles).
p2 <- plotGroups(
     ArchRProj = projSci, 
     groupBy = "Sample", 
     colorBy = "cellColData", 
     name = "TSSEnrichment",
     plotAs = "violin",
         alpha = 0.4,
         addBoxPlot = TRUE
    )
p2

# Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
     ArchRProj = projSci, 
     groupBy = "Sample", 
     colorBy = "cellColData", 
     name = "log10(nFrags)",
     plotAs = "ridges"
        )
p3
## Picking joint bandwidth of 0.05

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
p4

#To save editable vectorized versions of these plots, we use plotPDF().
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projSci, addDOC = FALSE, width = 4, height = 4)

#3.4 Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
#Because of how the data is stored and accessed, ArchR can compute fragment size distributions and TSS enrichment profiles from Arrow files very quickly.

#Fragment size distributions To plot the fragment size distributions of all samples, we use the plotFragmentSizes() function. Fragment size distributions in ATAC-seq can be quite variable across samples, cell types, and batches. Slight differences like those shown below are common and do not necessarily correlate with differences in data quality.
pfrag <- plotFragmentSizes(ArchRProj = projSci)
pfrag

#TSS enrichment profiles To plot TSS enrichment profiles, we use the plotTSSEnrichment() function. TSS enrichment profiles should show a clear peak in the center and a smaller shoulder peak right-of-center which is caused by the well-positioned +1 nucleosome.
pTSSEn <- plotTSSEnrichment(ArchRProj = projSci)
pTSSEn

#To save editable vectorized versions of these plots, we use plotPDF().
plotPDF(pfrag,pTSSEn, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projSci, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = projSci, outputDirectory = "SciTest", load = FALSE)




# Inferring Doublets
# After Arrow file creation, we can infer potential doublets (a single droplet
# containing multiple cells) that can confound downstream results. This is
# done using the addDoubletScores() function.

#commenting out doublet score becuase we saw this error message last time:
#ArchR logging to : ArchRLogs/ArchR-addDoubletScores-1f9f35c3dde9-Date-2021-07-22_Time-21-31-26.log
#If there is an issue, please report to github with logFile!
#2021-07-22 21:31:26 : Batch Execution w/ safelapply!, 0 mins elapsed.
#2021-07-22 21:31:26 : BAM_data (1 of 1) :  Computing Doublet Statistics, 0 mins elapsed.
#Warning: The following arguments are not used: row.names
#BAM_data (1 of 1) : UMAP Projection R^2 = 0.38719
#BAM_data (1 of 1) : Correlation of UMAP Projection is below 0.9 (normally this is ~0.99)
#This means there is little heterogeneity in your sample and thus doubletCalling is inaccurate.   
#force = FALSE, thus returning -1 doubletScores and doubletEnrichments!
#Set force = TRUE if you want to continue (not recommended).
#ArchR logging successful to : ArchRLogs/ArchR-addDoubletScores-1f9f35c3dde9-Date-2021-07-22_Time-21-31-26.log
#
#doubScores <- addDoubletScores(
#  input = ArrowFiles,
#  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
#  dimsToUse = 1:15 # Have to make upper dimension less than default of 30 since we only have 18 columns... 
#)

# Creating an ArchRProject 
# With our Arrow files in hand, we are now ready to create an ArchRProject. An ArchRProject is associated with a set of Arrow files and is the backbone of nearly all ArchR analyses.
#proj <- ArchRProject(
#  ArrowFiles = ArrowFiles, 
#  outputDirectory = "SciTest",
#  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
#)
### Using GeneAnnotation set by addArchRGenome(Hg19)!
### Using GeneAnnotation set by addArchRGenome(Hg19)!
### Validating Arrows…
### Getting SampleNames…
###
### Copying ArrowFiles to Ouptut Directory! If you want to save disk space set copyArrows = FALSE
### 1 2 3
### Getting Cell Metadata…
###
### Merging Cell Metadata…
### Initializing ArchRProject…
#
## We can also ask which data matrices are available within the ArchRProject which will be useful downstream once we start adding to this project:
#
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


##Now we can filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis.
#
#projSci <- filterDoublets(ArchRProj = projSci)
### Filtering 410 cells from ArchRProject!
### scATAC_BMMC_R1 : 243 of 4932 (4.9%)
### scATAC_CD34_BMMC_R1 : 107 of 3275 (3.3%)
### scATAC_PBMC_R1 : 60 of 2454 (2.4%)
#
## Dimensionality Reduction and Clustering
## ArchR implements an iterative LSI dimensionality reduction via the addIterativeLSI() function.
#
projSci <- addIterativeLSI(ArchRProj = projSci, useMatrix = "TileMatrix", name = "IterativeLSI")
#
## To call clusters in this reduced dimension sub-space, we use the addClusters() function which uses Seurat’s graph clustering as the default clustering method.
projSci <- addClusters(input = projSci, reducedDims = "IterativeLSI")
#
## Visualizing in a 2D UMAP Embedding
## We can visualize our scATAC-seq data using a 2-dimensional representation such as Uniform Manifold Approximation and Projection (UMAP). To do this, we add a UMAP embedding to our ArchRProject object with the addUMAP() function. This function uses the uwot package to perform UMAP.
#
cellColDataDF <- getCellColData(projSci)
write.csv(cellColDataDF, file='cell_column_data.csv')



projSci <- addUMAP(ArchRProj = projSci, reducedDims = "IterativeLSI")
### 16:32:38 UMAP embedding parameters a = 0.7669 b = 1.223
### 16:32:38 Read 10251 rows and found 30 numeric columns
### 16:32:38 Using Annoy for neighbor search, n_neighbors = 40
### 16:32:38 Building Annoy index with metric = cosine, n_trees = 50
### 0% 10 20 30 40 50 60 70 80 90 100%
### [—-|—-|—-|—-|—-|—-|—-|—-|—-|—-|
### **************************************************|
### 16:32:41 Writing NN index file to temp file /tmp/RtmpJ4Z9d9/file69ef223c61d2
### 16:32:41 Searching Annoy index using 10 threads, search_k = 4000
### 16:32:42 Annoy recall = 100%
### 16:32:44 Commencing smooth kNN distance calibration using 10 threads
### 16:32:45 Initializing from normalized Laplacian + noise
### 16:32:46 Commencing optimization for 200 epochs, with 623598 positive edges
### 16:32:58 Optimization finished
#
#projSciEmbeddingDF = getEmbedding(ArchRProj = projSci, embedding = "UMAP", returnDF = TRUE)
#write.csv(projSciEmbeddingDF, file='archr_umap_coords.csv')

projSciEmbeddingWClustersDF = getEmbedding(ArchRProj = projSci, embedding = "UMAP", returnDF = TRUE)
#message("row names")
#rownames(projSciEmbeddingWClustersDF)
#purchases$brand <- products$brand[products$id %in% purchases$product_id]
projSciEmbeddingWClustersDF$Clusters <- cellColDataDF$Clusters[match(row.names(cellColDataDF), row.names(projSciEmbeddingWClustersDF))]
write.csv(projSciEmbeddingWClustersDF, file='archr_umap_coords_clusters.csv')




## Using this UMAP, we can visualize various attributes of our cells which are stored in a matrix called cellColData in our ArchRProject. To do this, we use the plotEmbedding() function and we specify the variable to use for coloration via a combination of the colorBy and name parameters.
#
## For example, we can color by “Sample”:
#
pcellCollDataSampleUMAP <- plotEmbedding(ArchRProj = projSci, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
### ArchR logging to : ArchRLogs/ArchR-plotEmbedding-69ef25160449-Date-2020-04-21_Time-16-33-00.log
### If there is an issue, please report to github with logFile!
### Getting UMAP Embedding
### ColorBy = cellColData
### Plotting Embedding
### 1
### ArchR logging successful to : ArchRLogs/ArchR-plotEmbedding-69ef25160449-Date-2020-04-21_Time-16-33-00.log
#
## Or we can color by “Clusters”:
#
pCellCollDataClustersUMAP <- plotEmbedding(ArchRProj = projSci, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
### ArchR logging to : ArchRLogs/ArchR-plotEmbedding-69ef420f8c0a-Date-2020-04-21_Time-16-33-01.log
### If there is an issue, please report to github with logFile!
### Getting UMAP Embedding
### ColorBy = cellColData
### Plotting Embedding
### 1
### ArchR logging successful to : ArchRLogs/ArchR-plotEmbedding-69ef420f8c0a-Date-2020-04-21_Time-16-33-01.log
#
ggAlignPlots(pcellCollDataSampleUMAP, pCellCollDataClustersUMAP, type = "h")
#
#
## To save an editable vectorized version of this plot, we use the plotPDF() function.
#
plotPDF(pcellCollDataSampleUMAP, pCellCollDataClustersUMAP, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = projSci, addDOC = FALSE, width = 5, height = 5)
### [1] “plotting ggplot!”
### [1] “plotting ggplot!”
### [1] 0
#
## Download PDF : Plot-UMAP-Sample-Clusters.pdf
#
## Assigning Clusters with Gene Scores
## We can try to assign biological labels to these clusters using marker genes of known hematopoietic regulators. First, we add imputation weights using MAGIC to help smooth the dropout noise in our gene scores.
#
projSci <- addImputeWeights(projSci)
### ArchR logging to : ArchRLogs/ArchR-addImputeWeights-69ef433c71d0-Date-2020-04-21_Time-16-33-19.log
### If there is an issue, please report to github with logFile!
### 2020-04-21 16:33:19 : Computing Impute Weights Using Magic (Cell 2018), 0 mins elapsed.
#
## Now we can overlay our marker gene scores on our 2D UMAP embedding.
#
markerGenes  <- c(
    "CD34",  #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "MME", #B-Cell Trajectory
    "CD14", "MPO", #Monocytes
    "CD3D", "CD8A"#TCells
  )

pMarkerGenes <- plotEmbedding(
    ArchRProj = projSci, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(projSci)
)
## Getting ImputeWeights
## ArchR logging to : ArchRLogs/ArchR-plotEmbedding-69ef2f10a62d-Date-2020-04-21_Time-16-33-31.log
## If there is an issue, please report to github with logFile!
## Getting UMAP Embedding
## ColorBy = GeneScoreMatrix
## Getting Matrix Values…
## 2020-04-21 16:33:32 :
##
## Imputing Matrix
## Using weights on disk
## Using weights on disk
## Plotting Embedding
## 1 2 3 4 5 6 7 8 9
## ArchR logging successful to : ArchRLogs/ArchR-plotEmbedding-69ef2f10a62d-Date-2020-04-21_Time-16-33-31.log

# To plot a specific gene we can subset this plot list using the gene name.

p$CD14


# To plot all genes we can use cowplot to arrange the 9 different plots together. Each of these marker genes lights up the corresponding cell clusters. For example, we infer that the cells that have the highest gene score for CD3D, a known T cell marker, are in fact T cells.

#Rearrange for grid plotting
pMarkerGenesCow <- lapply(pMarkerGenes, function(x){
    x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),pMarkerGenesCow))

# To save an editable vectorized version of this plot, we use the plotPDF() function.

plotPDF(plotList = pMarkerGenes,
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
    ArchRProj = projSci,
    addDOC = FALSE, width = 5, height = 5)
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] “plotting ggplot!”
## [1] 0



# In addition to plotting gene scores per cell as a UMAP overlay, we can browse the local chromatin accessibility at these marker genes on a per cluster basis with genome browser tracks. To do this, we use the plotBrowserTrack() function which will create a list of plots, one for each of the genes specified by markerGenes.
pBrowserTrack <- plotBrowserTrack(
    ArchRProj = projSci,
    groupBy = "Clusters",
    geneSymbol = markerGenes,
    upstream = 50000,
    downstream = 50000
)

# To plot a track of a specific gene, we can simply select one from the list.
 grid::grid.newpage()
 grid::grid.draw(pBrowserTrack$CD14)




# We can save a multi-page PDF with a single page for each gene locus in our plot list using the plotPDF() function.

plotPDF(plotList = pBrowserTrack, 
    name = "Plot-Tracks-Marker-Genes.pdf", 
    ArchRProj = projSci, 
    addDOC = FALSE, width = 5, height = 5)
## NULL
## NULL
## NULL
## NULL
## NULL
## NULL
## NULL
## NULL
## NULL
## [1] 0

# Last but certainly not least, ArchR natively supports an interactive and dynamic genome browser that can be launched locally via a shiny app. To do this, we use the ArchRBrowser() function.

# ArchRBrowser(ArchRProj = projSci)
# This launches a dynamic genome browser session with a whole host of features including export of vectorized tracks for publication.


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
# First, lets remind ourselves of the cell types that we are working with in projHeme5 and their relative proportions.
table(projSci$Clusters)



# To identify marker genes based on gene scores, we call the getMarkerFeatures() function with useMatrix = "GeneScoreMatrix". We specify that we want to know the cluster-specific features with groupBy = "Clusters" which tells ArchR to use the “Clusters” column in cellColData to stratify cell groups.
markersGS <- getMarkerFeatures(
       ArchRProj = projSci, 
       useMatrix = "GeneScoreMatrix", 
       groupBy = "Clusters",
       bias = c("TSSEnrichment", "log10(nFrags)"),
       testMethod = "wilcoxon"
 )

markerGSList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerGSList$C1

# To visualize all of the marker features simultaneously, we can create a heatmap using the markerHeatmap() function, optionally supplying some marker genes to label on the heatmap via the labelMarkers parameter.
markersGSGenes  <- c(
      "CD34", #Early Progenitor
      "GATA1", #Erythroid
      "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
      "CD14", "CEBPB", "MPO", #Monocytes
      "IRF8", 
      "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )

heatmapGS <- markerHeatmap(
     seMarker = markersGS, 
     cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
     labelMarkers = markersGSGenes,
     transpose = TRUE
  )

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projSci, addDOC = FALSE)










# Often times, we are interested to know which peaks are unique to an individual cluster or a small group of clusters. 
# We can do this in an unsupervised fashion in ArchR using the addMarkerFeatures() function in combination with useMatrix = "PeakMatrix".
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

# We can use the getMarkers() function to retrieve particular slices of this SummarizedExperiment that we are interested in. The default behavior of this function is to return a list of DataFrame objects, one for each cell group.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

# If we are interested in the marker peaks for a specific cell group, we can access this from the list via the $ accessor.
#markerList$Erythroid

# Instead of a list of DataFrame objects, we can use getMarkers() to return a GRangesList object by setting returnGR = TRUE.
markers_gr <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
markers_gr
write.csv(markers_gr, file = "markers.csv")

# This GRangesList object can similarly be subset to a GRanges object for a particular cell group using the $ accessor.
markerList$C1


# ArchR provides multiple plotting functions to interact with the SummarizedExperiment objects returned by getMarkerFeatures().
# We can visualize these marker peaks (or any features output by getMarkerFeatures()) as a heatmap using the markerHeatmap() function.
#heatmapPeaks <- plotMarkerHeatmap(
#  seMarker = markersPeaks, 
#  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
#  transpose = TRUE
#  )
#
#
## We can plot this heatmap using draw().
##draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
#plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projSci, addDOC = FALSE)
#

# Instead of plotting a heatmap, we can also plot an MA or Volcano plot for any individual cell group. To do this, we use the plotMarkers() function. For an MA plot we specify plotAs = "MA". Here we specify the “Erythroid” cell group via the name parameter.
pma <- plotMarkers(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
# Similarly, for a Volcano plot, we specify plotAs = "Volcano".
pC1v <- plotMarkers(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pC1v
#
plotPDF(pma, pC1v, name = "C1-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projSci, addDOC = FALSE)
#
#
### Additionally we can see these peak regions overlayed on our browser tracks by passing the relevant peak regions to the features parameterin the plotBrowserTrack() function. This will add an additional BED-style track of marker peak regions to the bottom of our ArchR track plot. Here we specify plotting the GATA1 gene via the geneSymbol parameter
pPeakRegionsOverBrowserTrack <- plotBrowserTrack(
      ArchRProj = projSci,
      groupBy = "Clusters",
      geneSymbol = c("GATA1"),
      features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C1"],
      upstream = 50000,
      downstream = 50000
  )
#
##We can plot this using grid::grid.draw().
##grid::grid.newpage()
##grid::grid.draw(p$GATA1)
#i
plotPDF(pPeakRegionsOverBrowserTrack, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projSci, addDOC = FALSE)


# Here we perform a pairwise test between the “Erythroid” cell group and the “Progenitor” cell group.
markerTest <- getMarkerFeatures(
  ArchRProj = projSci, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1",
  bgdGroups = "C2"
)

#
# Motif Enrichment in Differential Peaks
projSci <- addMotifAnnotations(ArchRProj = projSci, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
       seMarker = markerTest,
       ArchRProj = projSci,
       peakAnnotation = "Motif",
       cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
    )

motifsUp

# To prepare this data for plotting with ggplot we can create a simplified data.frame object containing the motif names, the corrected p-values, and the significance rank.
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

# Using ggplot we can plot the rank-sorted TF motifs and color them by the significance of their enrichment. Here we use ggrepel to label each TF motif.
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
    ggrepel::geom_label_repel(
      data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + theme_ArchR() + 
   ylab("-log10(P-adj) Motif Enrichment") + 
   xlab("Rank Sorted TFs Enriched") +
   scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

# We can perform the same analyses for the peaks that are more accessible in the “Progenitor” cells by using peaks with Log2FC <= -0.5.
motifsDo <- peakAnnoEnrichment(
      seMarker = markerTest,
      ArchRProj = projSci,
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
    )

motifsDo

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
              data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
              size = 1.5,
              nudge_x = 2,
              color = "black"
        ) + theme_ArchR() + 
   ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
   scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "Erythroid-vs-Progenitor-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projSci, addDOC = FALSE)



# Saving and Loading an ArchRProject
# To easily save an ArchRProject for later use or for sharing with collaborators, we use the saveArchRProject() function. This copies the current ArchRProject object and all of the Arrow files to a specified directory. If we don’t specify an output directory (as below), saveArchRProject() uses the output directory that we specified upon creation of our ArchRProject. In this case that is the folder “HemeTutorial”.

projSci <- saveArchRProject(ArchRProj = projSci)
## Saving ArchRProject…
## Loading ArchRProject…
## Successfully loaded ArchRProject!

# When we are ready to load this saved ArchRProject we use the loadArchRProject() object and provide the path to the folder containing the saved ArchRProject object.

projSci <- loadArchRProject(path = "SciTest")
## Successfully loaded ArchRProject!

# Session Information
# This tutorial was run on the date specified below.

Sys.Date()
## [1] “2020-04-21”

# The sessionInfo() at run time was:
sessionInfo()

