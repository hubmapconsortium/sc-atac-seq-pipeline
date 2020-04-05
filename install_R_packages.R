install.packages(
	c(
		'devtools',
		'Matrix',
		'doSNOW',
		'plot3D',
		'optparse',
		'BiocManager'
	)
)
# Install SnapATAC
library(devtools)

install_github('r3fang/SnapATAC')
# Install R packages using BioConductor
# at https://bioconductor.org/install/#install-bioconductor-packages

# rtracklayer needed for cell selection work around described at
#      https://github.com/r3fang/SnapATAC/issues/139
# RFLPtools needed for the write.hclust call
# rmarkdown used when running an R markdown document
# dplyr needed for R markdown to create PDF files
BiocManager::install(
    c(
		'BSgenome.Hsapiens.UCSC.hg38',
		'rtracklayer',
		'RFLPtools',
		'rmarkdown',
		'dplyr',
		'motifmatchr',
		'SummarizedExperiment',
		'chromVAR'
	),
	ask=FALSE
)
