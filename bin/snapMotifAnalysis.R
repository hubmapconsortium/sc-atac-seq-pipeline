#!/usr/bin/env Rscript
library(SnapATAC)
library(optparse)
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.NCBI.GRCh38)

option_list = list(
  make_option(
    c("-s", "--snap_rds"),
    type="character",
    default=NULL,
    help="Snap object in RDS file"
  ),
  make_option(
    c("-p", "--snap_file"),
    type="character",
    default=NULL,
    help="Snap file"
  )
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$snap_rds)){
  print_help(opt_parser)
  stop("--snap_rds argument must be supplied (input file).n", call.=FALSE)
}

# Read RDS file to create the input snap object
x.sp = readRDS(opt$snap_rds)
# Add cell by peak matrix

# This code resets the path to the snap file to the current
# location. The original path was set when the snap object
# was created and may no longer be valid, especially when
# running under cwltool. addPmatToSnap uses this path to
# read the peaks matrix from the snap file, so it needs to
# be correct. The code to set the path is taken from
# createSnapSingle at
# https://rdrr.io/github/r3fang/SnapATAC/src/R/snap-utilities.R
x.sp@file <- rep(normalizePath(opt$snap_file), length(x.sp@barcode))

x.sp = addPmatToSnap(x.sp)

# Remove unwanted chromosomes, otherwise chromVAR throws the error with the message
# 'trying to load regions beyond the boundaries of non-circular sequence'
message(sprintf("Removing unwanted and mismatched chromosomes from peak matrix\n"))
# Provides seqlevels function
#chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM|chrUn", seqlevels(x.sp@feature))];

seqs_to_keep = as.vector(
  setdiff(
    intersect(x.sp@peak@seqnames, BSgenome.Hsapiens.NCBI.GRCh38@single_sequences@objnames),
    'MT'
  )
)

peaks_to_keep = x.sp@peak@seqnames %in% seqs_to_keep
x.sp = x.sp[, which(peaks_to_keep), mat="pmat"]

x.sp = makeBinary(x.sp, mat="pmat")

# SnapATAC also incorporates chromVAR (Schep et al) for motif variability analysis.

x.sp@mmat = runChromVAR(
    obj=x.sp,
    input.mat="pmat",
    genome=BSgenome.Hsapiens.NCBI.GRCh38,
    min.count=10,
    species="Homo sapiens"
)

message(sprintf("Writing the motif data to a Matrix Market format file\n"))

# Write cell by gene sparse matrix in Matrix Market format not CSV format
write.csv(
  x.sp@mmat,
  file="cellMotif.csv",
  na=''
)
#writeMM(obj = cellMotifData, file = "cellMotif.mtx")


#motif_i = "MA0497.1_MEF2C";
#dat = data.frame(x=x.sp@metaData[,"cluster"], y=x.sp@mmat[,motif_i]);
#p1 <- ggplot(dat, aes(x=x, y=y, fill=x)) +
#        theme_classic() +
#        geom_violin() +
#        xlab("cluster") +
#        ylab("motif enrichment") +
#        ggtitle(motif_i) +
#        theme(
#	      plot.margin = margin(5,1,5,1, "cm"),
#              axis.text.x = element_text(angle = 90, hjust = 1),
#              axis.ticks.x=element_blank(),
#	      legend.position = "none"
#   	      );
#
#plot_file_name=sprintf("Motif_%s_Enrichment_Per_Cluster", motif_i)
#ggsave(filename=plot_file_name, plot=p1)

