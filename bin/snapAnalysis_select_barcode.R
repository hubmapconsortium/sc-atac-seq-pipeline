#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-c", "--preferred_barcodes"), type="character", default=NULL,
              help="Genome mapping CSV file", metavar="character")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (!is.null(opt$preferred_barcodes)){
  #print(option_list);
  #print(opt$input_snap);

  #Step 1. Barcode selection
    #We select high-quality barcodes based on two criteria:
      #1) number of unique fragments
      #2) fragments in promoter ratio

  message(sprintf("Selecting barcodes from CSV file\n"));

  library(SnapATAC);

  barcodes = read.csv(
    opt$preferred_barcodes,
    head=TRUE
  );

  barcodes = barcodes[2:nrow(barcodes),];
  promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  UMI = log(barcodes$passed_filters+1, 10);
  data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
  barcodes$promoter_ratio = promoter_ratio;

  library(viridisLite);
  library(ggplot2);
  p1 = ggplot(
    data,
    aes(x= UMI, y= promoter_ratio)) +
    geom_point(size=0.1, col="grey") +
    theme_classic() +
    ggtitle("Barcode Selection UMI vs Promotor Ratio") +
    ylim(0, 1) + xlim(0, 6) +
    labs(x = "log10(UMI)", y="promoter ratio")

  ggsave("Barcode_Selection_UMI_vs_Promotor_Ratio.pdf", plot = p1);

  barcodes.sel = barcodes[which(UMI >= 3 & UMI <= 5 & promoter_ratio >= 0.15 & promoter_ratio <= 0.6),];
  rownames(barcodes.sel) = barcodes.sel$barcode;

  saveRDS(barcodes.sel, file = "selected_barcodes.rds")

} else {
  message(sprintf("No barcode CSV file provided\n"));
}

