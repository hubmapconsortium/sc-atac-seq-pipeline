#!/usr/bin/env Rscript
library(optparse)
library(SnapATAC)
library(rtracklayer)
library(GenomicRanges)
library(rhdf5)


ospj = function(...) {
  paste(c(...), collapse='/')
}

# TODO move/refactor this
supplementary_data_path = '/opt/supplementary-data'

default_gene_track = ospj(supplementary_data_path, 'gencode.v32.annotation.bed')
default_encode_blacklist = ospj(supplementary_data_path, 'hg38.blacklist.bed')
default_promoters = ospj(supplementary_data_path, 'hg38.promoters.bed')

snaptools_path = '/usr/local/bin/snaptools'
macs2_path = '/usr/local/bin/macs2'
# /TODO move/refactor this

option_list = list(
  make_option(
    c("-b", "--selected_barcodes"),
    type="character",
    default=NULL,
    help="Selected barcodes"
  ),
  make_option(
    c("-s", "--input_snap"),
    type="character",
    default=NULL,
    help="SNAP file name"
  ),
  make_option(
    c("-e", "--encode_blacklist"),
    type="character",
    default=default_encode_blacklist,
    help="ENCODE blacklist BED.GZ file"
  ),
  make_option(
    c("-g", "--gene_track"),
    type="character",
    default=default_gene_track,
    help="Gene track BED file"
  ),
  make_option(
    c("-a", "--gene_annotation"),
    type="character",
    default=NULL,
    help="Gene annotation GTF file"
  ),
  make_option(
    c("-p", "--promoters"),
    type="character",
    default=default_promoters,
    help="Promoters BED file"
  ),
  make_option(
    c("-n", "--processes"),
    type="integer",
    default=1,
    help="Number of subprocesses/threads to use"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$input_snap)){
  print_help(opt_parser)
  stop("--input_snap argument must be supplied (input file).", call.=FALSE)
}

mytmpdir <- Sys.getenv("TMP")
sprintf("TMP is: %s", mytmpdir)
mytmpdir <- Sys.getenv("TMPDIR")
sprintf("TMPDIR is: %s", mytmpdir)
mytmpdir <- Sys.getenv("TEMP")
sprintf("TEMP is: %s", mytmpdir)
mytmpdir <- tempdir()
sprintf("tempdir() is: %s", mytmpdir)

#write("TMPDIR = D:/mnt/", file=file.path(Sys.getenv('TMPDIR'), '.Renviron'))

x.sp = createSnap(
  file=opt$input_snap,
  sample="Input_snap",
  num.cores=opt$processes
)
x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=opt$processes)
# Assign column (bin) names in the sparse cell-by-bin matrix
dimnames(x.sp@bmat)[[2]] = x.sp@feature$name

if (!is.null(opt$selected_barcodes)) {
  message(sprintf("Reading selected barcode file\n"))
  barcodes.sel <- readRDS(opt$selected_barcodes)
  x.sp = x.sp[which(x.sp@barcode %in% barcodes.sel$barcode),]
  x.sp@metaData = barcodes.sel[x.sp@barcode,]

} else if (!is.null(opt$gene_annotation)) {
  # If there is no barcode CSV file
  # Otherwise select cells based on "coverage" and "promoter ratio -"
  # percentage of reads mapped to promoters per cell
  # https://github.com/r3fang/SnapATAC/issues/139

  message("Selecting barcodes base on coverage and promoter ratio")
  gtf.gr <- rtracklayer::import(opt$gene_annotation)
  gene.gr <- gtf.gr[gtf.gr$type == "gene"]

  # extract promoter region for each gene
  promoter.gr <- reduce(promoters(gene.gr, upstream=2000, downstream=0))

  ov = findOverlaps(x.sp@feature, promoter.gr)

  # find promoter overlapping bins
  idy = queryHits(ov)
  log_cov = log10(SnapATAC::rowSums(x.sp, mat="bmat")+1)
  promoter_ratio = Matrix::rowSums(x.sp@bmat[,idy]) / Matrix::rowSums(x.sp@bmat)

  pdf(file="PromotorRatioLogPlot.pdf")
  plot(log_cov, promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio (promotor)", ylim=c(0,1 ))
  dev.off()

  # From Matt: Don't select cutoff 2/24/2020
  # Choose the cutoff based on the plot
  # however we need column idx for ChromVAR motif analysis...
  #idx = which(promoter_ratio > 0.2 & promoter_ratio < 0.8 & log_cov > 3)
  #x.sp = x.sp[idx,]

} else {

  if (is.null(opt$promoters)){
    print_help(opt_parser)
    stop("--promotors argument must be supplied if no barcode CSV or gene annotation file is provided", call.=FALSE)
  }

  # From Dihn's KC20_Test.Rmd file
  # The input for example is "hg38.promoters.bed"
  calBmatCor(x.sp)

  promoter.df = read.table(opt$promoters)
  promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]))
  ov = findOverlaps(x.sp@feature, promoter.gr)
  idy = queryHits(ov)
  promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat")

  pdf(file="PromotorRatioLogPlot.pdf")
  plot(
    x=log(SnapATAC::rowSums(x.sp, mat="bmat") + 1,10),
    y=promoter_ratio,
    cex=0.5,
    col="grey",
    xlab="log(count)",
    ylab="FIP Ratio (promoter)",
    ylim=c(0, 1)
  )
  dev.off()

  # TODO write the promoter ratio to a column in another output csv or its own csv
}

message(sprintf("Writing bin data to csv\n"))
Bins <- x.sp@feature
write.csv(Bins, file = "Bins.csv", row.names = FALSE)

# Filter cells based on attributes
# https://rdrr.io/github/r3fang/SnapATAC/man/filterCells.html
plotBarcode(x.sp, pdf.file.name = "BarcodeQualityControlDistributionBefore.pdf", col="grey", border="grey")
x.sp = filterCells(
    obj=x.sp,
    subset.names=c("fragment.num", "UMI", "mito.ratio"),
    low.thresholds=c(2000, 1000, 0),
    high.thresholds=c(Inf, Inf, 0.5)
)
plotBarcode(x.sp, pdf.file.name = "BarcodeQualityControlDistributionAfter.pdf", col="grey", border="grey")

calBmatCor(x.sp)

# We will convert the cell-by-bin count matrix to a binary matrix. Some items
# in the count matrix have abnormally high coverage perhaps due to the alignment
# errors. Therefore, we next remove 0.1% items of the highest coverage in the
# count matrix and then convert the remaining non-zero items to 1.
message(sprintf("Converting cell-by-bin-matrix to a binary matrix\n"))
x.sp = makeBinary(x.sp, mat="bmat")

# The bin coverage roughly obeys a log normal distribution. We remove the
# top 5% bins that overlap with invariant features such as promoters of the house keeping genes.
message(sprintf("Removing top 5 percent of bins that overlap with invariant features\n"))
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)

# Filter out bins, i.e. bins with less than 5% or more than 95% coverage
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
x.sp = x.sp[, idy, mat="bmat"]


message("Creating coverage histogram")
pdf(file="CoverageHistogram.pdf")
hist(
  bin.cov[bin.cov > 0],
  xlab="log10(bin cov)",
  main="log10(Bin Cov)",
  col="lightblue",
  xlim=c(0, 5)
)
dev.off()

# Bin filtering
# Filter out any bins overlapping with the ENCODE blacklist to prevent from potential artifacts.
message("Doing bin filtering")

# Read the black list table and ignore missing columns
# NOTE: ignoring missing columns could cause a problem
# if by some chance the chr region start or end was missing from the BED file
black_list = read.table(opt$encode_blacklist, header=FALSE, row.names=NULL, fill = TRUE)

black_list.gr = GRanges(
  black_list[,1],
  IRanges(black_list[,2], black_list[,3])
)

#idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))
#if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}

# Second, we remove unwanted chromosomes.
message("Removing unwanted random and M chromosomes")
chrrandomM.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))]
chrrandomM.exclude
idy = grep(paste(chrrandomM.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}

#message(sprintf("Removing unwanted decoy chromosomes\n"))
#chrdecoy.exclude = seqlevels(x.sp@feature)[grep("decoy", seqlevels(x.sp@feature))]
#chrdecoy.exclude
#idydecoy = grep(paste(chrdecoy.exclude, collapse="|"), x.sp@feature)
#if(length(idydecoy) > 0){x.sp = x.sp[,-idydecoy, mat="bmat"]}
#
#message(sprintf("Removing unwanted Unknown chromosomes\n"))
#chrUn.exclude = seqlevels(x.sp@feature)[grep("chrUn", seqlevels(x.sp@feature))]
#chrUn.exclude
#idyUnknown = grep(paste(chrUn.exclude, collapse="|"), x.sp@feature)
#if(length(idyUnknown) > 0){x.sp = x.sp[,-idyUnknown, mat="bmat"]}

# from KC20_Test.Rmd Dinh's script
# however unwanted chromosomes and black listed items are not actually
# removed in her script.
# TODO: Last word is that they should be removed from the BAM file, not sure about unwanted chromosomes
#idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr))
#idy2 = grep("chrM|random", x.sp@feature)
#idy = unique(c(idy1, idy2))


# Dimensionality reduction
# We compute diffusion maps for dimentionality reduction.

# Remove empty rows from bmat
# so runDiffusionMaps does not throw an error
# https://stackoverflow.com/questions/6437164/removing-empty-rows-of-a-data-file-in-r
idempty <- which(Matrix::rowSums(x.sp@bmat) == 0)

if(length(idempty) > 0){
    message("Removing rows with all zeros from bmat")
    x.sp <- x.sp[-idempty,, mat='bmat']
}

write.table(dimnames(x.sp@bmat)[[1]], 'barcodes.txt', col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(dimnames(x.sp@bmat)[[2]], 'bins.txt', col.names=FALSE, row.names=FALSE, quote=FALSE)
writeMM(x.sp@bmat, 'filtered_cell_by_bin.mtx')

message("Computing diffusion maps")
x.sp = runDiffusionMaps(
  obj=x.sp,
  input.mat="bmat",
  num.eigs=50
)

# This step is from https://github.com/r3fang/SnapATAC/blob/master/examples/10X_snATAC/README.md
#row.covs = log10(Matrix::rowSums(x.sp@bmat)+1)
#message(sprintf("Executing step 5: row.covs:\n"))
#head(row.covs)
#
#row.covs.dens = density(
#      x = row.covs,
#      bw = 'nrd', adjust = 1
#)
#
#message(sprintf("Executing step 5: row.covs.dens:\n"))
#head(row.covs.dens)
#
#
#sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
#set.seed(1)
#
#message(sprintf("Executing step 5: nrow(x.sp):\n"))
#nrow(x.sp)
#
##idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 10000, prob = sampling_prob))
#idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = round(nrow(x.sp)*.75), prob = sampling_prob))
#x.landmark.sp = x.sp[idx.landmark.ds,]
#x.query.sp = x.sp[-idx.landmark.ds,]
#
#x.landmark.sp = runDiffusionMaps(
#       obj= x.landmark.sp,
#       input.mat="bmat",
#       num.eigs=50
#     )
#
#x.query.sp = runDiffusionMapsExtension(
#     obj1=x.landmark.sp,
#     obj2=x.query.sp,
#     input.mat="bmat"
#   )
#
#x.landmark.sp@metaData$landmark = 1
#x.query.sp@metaData$landmark = 0
#x.sp = snapRbind(x.landmark.sp, x.query.sp)
### combine landmarks and query cells
#x.sp = x.sp[order(x.sp@sample),]; # IMPORTANT
#rm(x.landmark.sp, x.query.sp); # free memory
#

# Determine significant components
# We next determine the number of reduced dimensions to include for downstream
# analysis. We use an ad hoc method by simply looking at a pairwise plot and
# select the number of dimensions in which the scatter plot starts looking like
#a blob. In the below example, we choose the first 20 dimensions.
#message(sprintf("Creating eigen plots\n"))
plotDimReductPW(
  obj=x.sp,
  eigs.dims=1:50,
  point.size=0.3,
  point.color="grey",
  point.shape=19,
  point.alpha=0.6,
  down.sample=5000,
  pdf.file.name='EigenPlots.pdf',
  pdf.height=7,
  pdf.width=7
)

# Graph-based clustering
# Using the selected significant dimensions, we next construct a K Nearest Neighbor
# (KNN) Graph (k=15). Each cell is a node and the k-nearest neighbors of each cell
# are identified according to the Euclidian distance and edges are draw between
# neighbors in the graph.
message(sprintf("Graph-based clustering\n"))
x.sp = runKNN(
  obj=x.sp,
  eigs.dims=1:20,
  k=15
)
x.sp=runCluster(
  obj=x.sp,
  tmp.folder=tempdir(),
  louvain.lib="R-igraph",
  seed.use=10
)
x.sp@metaData$cluster = x.sp@cluster

cellClusterAssignment <- as.data.frame(x.sp@cluster)
# Add a barcode column
cellClusterAssignment$BarcodeID <- dimnames(x.sp@bmat)[[1]]

# Update: only including the numeric barcode IDs in most files, since most people won't
# care about the actual barcodes for most uses, and including these in every single
# file just has the effect of inflating the file sizes. If people actually do care about
# the barcode sequence of a cell, they'll still be able to use 'cellBarcodes.csv' to obtain it.
# Look up the barcode string in the look up table using the barcode ID
# column and add a new column for the barcode string
#cellClusterAssignment$Barcode <- BarcodeLookupTable[match(cellClusterAssignment$BarcodeID, BarcodeLookupTable$BarcodeID), "Barcode"]
names(cellClusterAssignment)[1] <- 'Cluster'
# write it to csv and only quote the Barcode in column 3
#write.csv(cellClusterAssignment, file = "cellClusterAssignment.csv", quote = c(3), row.names = FALSE)
message("Writing cell cluster assignment csv")
write.csv(cellClusterAssignment, file = "cellClusterAssignment.csv", quote = FALSE, row.names = FALSE)

## Step 8. Visualization
message("Creating embedding UMAP plot")
x.sp = runViz(
  obj=x.sp,
  tmp.folder=tempdir(),
  dims=2,
  eigs.dims=1:20,
  method="umap",
  seed.use=10
)

#par(mfrow = c(2, 2))
plotViz(
  pdf.file.name='Embedding_UMAP.pdf',
  obj=x.sp,
  method="umap",
  main="Embedding UMAP",
  point.color=x.sp@cluster,
  point.size=1,
  point.shape=19,
  point.alpha=0.8,
  text.add=TRUE,
  text.size=1.5,
  text.color="black",
  text.halo.add=TRUE,
  text.halo.color="white",
  text.halo.width=0.2,
  down.sample=10000,
  legend.add=FALSE
)

message("Saving UMAP coordinates and cluster assignments as CSV")
dim_reduced_data = data.frame(x.sp@umap, row.names=x.sp@barcode)
dim_reduced_data$cluster = as.numeric(x.sp@cluster)
write.csv(dim_reduced_data, file='umap_coords_clusters.csv')

message(paste("Reading gene track from", opt$gene_track))
genes = read.table(opt$gene_track)
genes.gr = GRanges(genes[,1],
  IRanges(genes[,2], genes[,3]), name=genes[,4]
)

message("Writing gene ranges to csv")
genes_df = as(genes.gr, "data.frame")
write.csv(genes_df, file = "GenesRanges.csv", row.names = FALSE)

message("Creating cell by gene matrix")
x.sp = createGmatFromMat(
  obj=x.sp,
  input.mat="bmat",
  genes=genes.gr,
  do.par=TRUE,
  num.cores=opt$processes
)

raw_cell_by_gene_filename = 'cell_by_gene_raw.mtx'
message(paste("Writing raw cell by gene data to", raw_cell_by_gene_filename))
writeMM(x.sp@gmat, raw_cell_by_gene_filename)

# smooth the cell-by-gene matrix
message("Smoothing cell-by-gene matrix")
x.sp = runMagic(
  obj=x.sp,
  input.mat="gmat",
  step.size=3
)

smooth_cell_by_gene_filename = 'cell_by_gene_smoothed.hdf5'
message(paste("Writing smoothed cell by gene data to", smooth_cell_by_gene_filename))
h5createFile(smooth_cell_by_gene_filename)
h5write(as.matrix(x.sp@gmat), smooth_cell_by_gene_filename, 'cell_by_gene_smoothed', level=0)
h5write(dimnames(x.sp@gmat)[[2]], smooth_cell_by_gene_filename, 'genes')
h5write(x.sp@barcode, smooth_cell_by_gene_filename, 'barcodes')

# Step 11. Identify peaks
# Next we aggregate cells from the each cluster to create an ensemble track for
# peak calling and visualization. This step will generate a .narrowPeak file that
# contains the identified peak and .bedGraph file for visualization. To obtain
# the most robust result, we don't recommend to perform this step for clusters
# with cell number less than 100. In the below example, SnapATAC creates
# atac_v1_adult_brain_fresh_5k.1_peaks.narrowPeak and
# atac_v1_adult_brain_fresh_5k.1_treat_pileup.bdg. bdg file can be compressed
# to bigWig file using bedGraphToBigWig for IGV or Genome Browser visulization.
message("Identifying peaks")

# Identify peaks using selected cells. Fragments belonging to a subset of cells
# are extracted and used to identify peaks using MACS2. This function requires
# "MACS2" and "snaptools" preinstalled and excutable.
# https://rdrr.io/github/r3fang/SnapATAC/man/runMACS.html
runMACS(
  obj=x.sp,
  output.prefix="overall",
  path.to.snaptools=snaptools_path,
  path.to.macs=macs2_path,
  #gsize: effective genome size. 'hs' for human, 'mm' for mouse, 'ce' for C. elegans, 'dm' for fruitfly (default: None)
  gsize="hs",
  buffer.size=500,
  num.cores=opt$processes,
  macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
  tmp.folder=tempdir()
)

# call peaks for all cluster with more than 100 cells
cluster_counts = table(x.sp@cluster)
clusters_sel = names(cluster_counts)[which(cluster_counts > 100)]
peaks.ls = lapply(
  clusters_sel,
  function(cluster) {
    message(paste('Running MACS2 for cluster', cluster))
    runMACS(
      obj=x.sp[which(x.sp@cluster==cluster),],
      output.prefix=paste0('cluster_', cluster),
      path.to.snaptools=snaptools_path,
      path.to.macs=macs2_path,
      gsize="hs", # mm, hs, etc
      buffer.size=500,
      num.cores=opt$processes,
      macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
      tmp.folder=tempdir()
    )
  }
)

# TODO: possibly construct file names ourselves since we know exactly which
# ones should exist, given the 'clusters_sel' vector above
cluster_peak_files = list.files(pattern='cluster_\\d+_peaks\\.narrowPeak')
if (length(clusters_sel) == 0L) {
    message("No clusters with enough cells; using overall peaks")
    cluster_peak_files = list.files(pattern='overall_peaks.narrowPeak')
}

message("Cluster peak files:")
cluster_peak_files

peak.gr.ls = lapply(
  cluster_peak_files,
  function(filename) {
    peak.df = read.table(filename, col.names=as.character(1:10))
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  }
)
peak.gr = reduce(Reduce(c, peak.gr.ls))

message("Creating a cell by peak matrix")

peaks.df = as.data.frame(peak.gr)[,1:3]

write.table(
  peaks.df,
  file="peaks.combined.bed",
  quote=FALSE,
  sep="\t",
  dec=".",
  row.names=FALSE,
  col.names=FALSE
)

write.csv(peaks.df, file = "peaksAllCells.csv", row.names = FALSE)

saveRDS(x.sp, file="peaks_snap.rds")
