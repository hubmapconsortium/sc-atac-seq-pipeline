# Overview

The HuBMAP Consortium sc-atac-seq pipeline is a three-stage pipeline for
analyzing scATAC-seq data sets, composed of [SnapTools](https://github.com/r3fang/SnapTools "SnapTools"),
[SnapATAC](https://github.com/r3fang/SnapATAC "SnapATAC"),
and [chromVAR](https://bioconductor.org/packages/release/bioc/html/chromVAR.html "chromVAR").
Source code can be found at https://github.com/hubmapconsortium/sc-atac-seq-pipeline

SnapTools performs quantification using a specified aligner,
and HuBMAP has standardized on BWA with the GRCh38 reference genome.
SnapTools divides the genome into non-overlapping bins of user-specified size
(we use 1000, 5000, and 10000), produces FASTQC analysis of the input fastq files,
and produces a binary cell-by-bin matrix
denoting whether reads in each cell were aligned to each bin.

The SnapATAC
secondary analysis pipeline filters bins based on the ENCODE genomic region
blacklist, performs t-SNE dimensionality reduction, and selects peaks from
all available bins. The chromVAR tool performs motif analysis, assigns motifs
to transcription factors, and computes differential enrichment of transcription
factors across cells in the data set.

# Requirements
Running the pipeline requires a CWL workflow execution engine, and we recommend
the cwltool reference implementation, which is written in Python. This can be
installed in a sufficiently recent Python environment with pip install
[cwltool](https://github.com/common-workflow-language/cwltool),
after which the pipeline can be invoked as:

`cwltool create_snap_and_analyze.cwl create_snap_and_analyze.json`

# Supplementary Data
The HuBMAP sc-atac-seq pipeline uses the Genome Reference Consortium human genome,
build 38 (GRCh38).  A BWA generated set of index files is required for the
reference genome. Using an alternate reference or index is not currently
supported without rebuilding the sc-atac-seq Docker container, though one can
build an alternate container by modifying the Dockerfile.

# Inputs
## Required
* sequence_directory\
A directory for the pipeline to search for fastq or fastq.gz files. The pipeline
only works on paired end reads and expects, for historical reasons, the paired
end read files to be named
`<some_name>*_R1*.fastq` and `<some_name>*_R3*.fastq`. If a file
containing barcodes `<some_name>*_R2*.fastq` is found the barcodes will be
read and added to the read IDs in the paired end fastq files

* input_reference_genome\
A fasta file of the GRCh38 reference genome

* encode_blacklist\
The genome BED file with blacklisted regions. If reads overlap with blacklisted
regions they are removed from the BAM file. E.g `hg38.blacklist.bed`

* gene_track\
A BED file of gene tracks. Used to create the cell by gene matrix.
E.g `gencode.v32.annotation.bed`

## Optional

* reference_genome_index\
A .gz file containing the BWA generated index of the GRCh38 reference genome.
I.e. the ".bwt", ".sa", ".ann", ".pac", ".amb" files generated by BWA
indexing. If this file is provided the index will not have to be generated
by the pipeline saving some time.

* preferred_barcodes\
A CSV file. Cells will be selected based on the barcodes listed in the CSV file. E.g
[barcodes CSV](http://renlab.sdsc.edu/r3fang/share/github/Mouse_Brain_10X/atac_v1_adult_brain_fresh_5k_singlecell.csv
 "barcode CSV")

* gene_annotation\
A GTF file of gene annotations. If there is no preferred_barcodes CSV file
 cells will be selected based on "coverage" and "promoter ratio -"
i.e. the percentage of reads mapped to promoters per cell

* promoters\
A BED file of promoter regions used to calculate a promoter ratio. If no
preferred barcode or gene annotation file is provided this option is required.


# Outputs

* Bins.csv\
A CSV file providing sequence name and bin information

* cellBarcodes.CSV\
A CSV file with barcode ID and barcode

* cellByBin_summary.csv\
A CSV file with barcode ID and bin number

* cellClusterAssignment.csv\
A CSV file with barcode ID and cluster number

* GenesRanges.csv\
A CSV file providing sequence, gene name and gene location information

* cellByGene.mtx\
A file with the cell by gene matrix in Matrix Market format

* cellGenes.csv\
A CSV file with gene ID and gene name

* peaksAllCells.csv\
A CSV file with sequence name and peak start and end
