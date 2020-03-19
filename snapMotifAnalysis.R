#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-s", "--snap_file"), type="character", default=NULL,
              help="Selected barcodes", metavar="character")

  # ,make_option(c("-r", "--tmpdir"), type="character", default=NULL,
  #            help="Temporary directory path", metavar="character")

  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$snap_file)){
  print_help(opt_parser)
  stop("--snap_file argument must be supplied (input file).n", call.=FALSE)
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

library(SnapATAC);
x.sp = createSnap(
  file=opt$snap_file,
  sample="snap file with peaks data",
  num.cores=6
);


# Add cell by peak matrix
x.sp = addPmatToSnap(x.sp);
x.sp = makeBinary(x.sp, mat="pmat");


# SnapATAC also incorporates chromVAR (Schep et al) for motif variability analysis.
library(chromVAR);
library(motifmatchr);
library(SummarizedExperiment);
library(BSgenome.Hsapiens.UCSC.hg38);


x.sp@mmat = runChromVAR(
    obj=x.sp,
    input.mat="pmat",
    genome=BSgenome.Hsapiens.UCSC.hg38,
    min.count=10,
    species="Homo sapiens"
);

message(sprintf("Writing the cell by gene data to a Matrix Market format file\n"))
cellMotifData <- x.sp@mmat;
cellMotifMatrix <- as.matrix(cellMotifData)
cellMotifFrame <- as.data.frame(cellMotifMatrix);

# Write cell by gene sparse matrix in Matrix Market format not CSV format
write.csv(cellMotifFrame, file = "cellMotif.csv");
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
#
#
#motif_i = "MA0660.1_MEF2B";
#dat = data.frame(x=x.sp@metaData[,"cluster"], y=x.sp@mmat[,motif_i]);
#p2 <- ggplot(dat, aes(x=x, y=y, fill=x)) + 
#	  theme_classic() +
#          geom_violin() + 
#          xlab("cluster") +
#          ylab("motif enrichment") + 
#          ggtitle(motif_i) +
#          theme(
#              plot.margin = margin(5,1,5,1, "cm"),
#              axis.text.x = element_text(angle = 90, hjust = 1),
#              axis.ticks.x=element_blank(),
#	      legend.position = "none"
#          );
#
#plot_file_name=sprintf("Motif_%s_Enrichment_Per_Cluster", motif_i)
#ggsave(filename=plot_file_name, plot=p2)

