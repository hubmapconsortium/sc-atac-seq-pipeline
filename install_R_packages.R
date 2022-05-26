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

# Docker image rocker/tidyverse includes devtools so no need to install it
# First, install devtools (for installing GitHub packages) if it isn’t already installed:
#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.13", ask=FALSE)

tryCatch({
    BiocManager::install("DirichletMultinomial", ask=FALSE)
},
    error = function(e) {
    message("Error installing DirichletMultinomial")
    message(e$message)
    quit("no", -1)
  }
)

# Then, install ArchR:
tryCatch({
    devtools::install_github("GreenleafLab/ArchR@b9ee2663d9ba7d58c6737a7a8bf2b3614bf26866", repos = BiocManager::repositories())
},
    error = function(e) {
    message("Error installing GreenleafLab/ArchR")
    message(e$message)
    quit("no", -1)
  }
)

## Lastly, install all of the ArchR dependencies that arent installed by default:
library(ArchR)

tryCatch({
    ArchR::installExtraPackages()
},
    error = function(e) {
    message("Error installing ArchR extra packages")
    message(e$message)
    quit("no", -1)
  }
)

#
## If any of these steps fails, you should identify the offending package and
## troubleshoot that individual installation before proceeding. The one exception
## is Cairo (see below) which is installed by the ArchR::installExtraPackages()
## function. Cairo is not required but is highly recommended.

# We installed MACS2 using python and requirements.txt
## It is also highly recommended that you install MACS2, which requires python,
## and have the macs2 executable in your PATH variable. This will allow ArchR to call peaks using MACS2.
##pip install MACS2
#
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
