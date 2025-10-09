install.packages(
c(
  'optparse'
  ),
 Ncpus=6
)

# ArchR is designed to be run on Unix-based operating systems such as macOS and linux.
# ArchR is NOT supported on Windows or other operating systems.

# ArchR installation currently requires devtools and BiocManager for installation
# of GitHub and Bioconductor packages. Run the following commands to install the
# various dependencies used by ArchR:

# Docker image greenleaflab/archr:1.0.3-base-r4.4 includes devtools, BiocManager, and ArchR so no need to install them

# If any of these steps fails, you should identify the offending package and
# troubleshoot that individual installation before proceeding. The one exception
# is Cairo (see below) which is installed by the ArchR::installExtraPackages()
# function. Cairo is not required but is highly recommended.

# We installed MACS2 using python and requirements.txt

tryCatch({
    devtools::install_github("GreenleafLab/chromVARmotifs")
},
    error = function(e) {
    message("Error installing GreenleafLab chromVARmotifs")
    message(e$message)
    quit("no", -1)
  }
)

tryCatch({
    install.packages("magick", version = "2.7.3")
},
    error = function(e) {
    message("Error installing magick")
    message(e$message)
    quit("no", -1)
  }
)

tryCatch({
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", ask=FALSE)
},
    error = function(e) {
    message("Error installing UCSC hg38")
    message(e$message)
    quit("no", -1)
  }
)
