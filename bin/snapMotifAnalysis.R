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
message(sprintf("Removing unwanted chromosomes from PMAT\n"))
# Provides seqlevels function
#chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM|chrUn", seqlevels(x.sp@feature))];
chr.exclude = seqlevels(x.sp@feature)[grep("random|MT", seqlevels(x.sp@feature))]
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="pmat"]}

message(sprintf("Making PMAT binary\n"))
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
cellMotifData <- x.sp@mmat
cellMotifMatrix <- as.matrix(cellMotifData)
cellMotifFrame <- as.data.frame(cellMotifMatrix)

# Write cell by gene sparse matrix in Matrix Market format not CSV format
write.csv(cellMotifFrame, file = "cellMotif.csv")
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

