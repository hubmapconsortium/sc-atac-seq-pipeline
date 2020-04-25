library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
peak_file = args[1]
bam_list = args[2:length(args)]

read_Narrowpeaks = function (filename, width=500, non_overlapping = TRUE)
{
  cn <- c("chr", "start", "end", "name", "score", "strand",
          "fc", "pval", "qval", "summit")
  bed <- read.delim(file = filename, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE, col.names = cn)
  bed[, "summit"] <- bed[, "start"] + bed[, "summit"]
  bed <- as(bed, "DataFrame")
  bed <- makeGRangesFromDataFrame(bed[, c("chr", "summit",
                                          "score", "qval", "name")], start.field = "summit", end.field = "summit",
                                  keep.extra.columns = TRUE)
  bed <- resize(bed, width = width, fix = "center")
  if (non_overlapping) {
    bed <- sortSeqlevels(bed)
    bed <- sort(bed)
    keep_peaks <- seq_along(bed)
    while (!(isDisjoint(bed[keep_peaks]))) {
      chr_names <- as.character(seqnames(bed[keep_peaks]))
      starts <- start(bed[keep_peaks])
      ends <- end(bed[keep_peaks])
      overlap_next <- intersect(
        which(
          chr_names[seq_len(length(keep_peaks) - 1)]
          == chr_names[seq_len(length(keep_peaks) - 1) + 1]
        ),
        which(
          ends[seq_len(length(keep_peaks) - 1)]
          >= starts[seq_len(length(keep_peaks) - 1) + 1]
        )
      )
      overlap_previous <- overlap_next + 1
      overlap_comparison <- bed[keep_peaks[overlap_previous]]$qval >
        bed[keep_peaks[overlap_next]]$qval
      discard <- keep_peaks[c(overlap_previous[!overlap_comparison],
                              overlap_next[overlap_comparison])]
      keep_peaks <- keep_peaks[!(keep_peaks %in% discard)]
    }
    bed <- bed[keep_peaks]
  }
  return(bed)
}

cat('Reading peaks\n')
peaks = read_Narrowpeaks(peak_file)
peaks = peaks[seqnames(peaks) %in% seqnames(BSgenome.Hsapiens.NCBI.GRCh38)]

samples = sapply(strsplit(bam_list, '_'), function(x) {x[1]})

cat(paste(c('Reading ', length(bam_list), ' BAM files\n'), collapse=''))
fragment_counts = getCounts(
  bam_list,
  peaks,
  paired=TRUE,
  by_rg=FALSE,
  format="bam",
  colData = DataFrame(sample=samples)
)
cat('Computing GC bias\n')
fragment_counts <- addGCBias(fragment_counts, genome=BSgenome.Hsapiens.NCBI.GRCh38)

cat('Loading motifs\n')
motifs <- getJasparMotifs()
cat('Matching motifs\n')
counts_filtered = filterPeaks(fragment_counts)
motif_ix <- matchMotifs(motifs, counts_filtered, genome=BSgenome.Hsapiens.NCBI.GRCh38)
cat(typeof(motif_ix))
#If possible, output this as well

cat('Computing deviations\n')
dev <- computeDeviations(object=counts_filtered, annotations=motif_ix)

scores = deviationScores(dev)
cat(typeof(scores))
#This should be a numerical matrix.  Output it somehow
write.csv(scores)
