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
message("Removing unwanted and mismatched chromosomes from peak matrix")
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

rse <- SummarizedExperiment(
  assays = list(counts = t(x.sp@pmat)),
  rowRanges = x.sp@peak,
  colData = DataFrame(Cell_Type=1:nrow(x.sp@pmat), depth=Matrix::rowSums(x.sp@pmat))
)
rse <- addGCBias(rse, genome=BSgenome.Hsapiens.NCBI.GRCh38)
motifs <- getJasparMotifs(collection = "CORE")
motif_mm <- matchMotifs(motifs, rse, genome=BSgenome.Hsapiens.NCBI.GRCh38)
dev <- computeDeviations(object=rse, annotations=motif_mm)
dev_mat = t(assay(dev))

deviation_scores = deviationScores(dev)
deviation_score_filename = 'chromvar_deviation_scores.csv'
message(paste('Saving deviation scores to', deviation_score_filename))
write.csv(deviation_scores, file=deviation_score_filename, na='')

variability_scores = computeVariability(dev)
variability_score_filename = 'chromvar_variability_scores.csv'
message(paste('Saving deviation scores to', variability_score_filename))
write.csv(variability_scores, file=variability_score_filename, na='')

pdf('output.pdf', width=7, height=16)
heatmap(scores, col=brewer.pal(11, 'RdBu'), scale='none', na.rm=TRUE)
dev.off()

save(
  rse,
  motif_mm,
  dev,
  dev_mat,
  deviation_scores,
  variability_scores,
  file='chromvar_data.RData'
)

message("Writing the motif data to CSV")

write.csv(
  dev_mat,
  file="cellMotif.csv",
  na=''
)
