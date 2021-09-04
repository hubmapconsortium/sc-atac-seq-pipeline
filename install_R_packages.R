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

tryCatch({
    BiocManager::install("DirichletMultinomial")
},
    error = function(e) {
    message("Error installing DirichletMultinomial")
    message(e$message)
    quit("no", -1)
  }
)

# Then, install ArchR:
#devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
tryCatch({
    devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
},
    error = function(e) {
    message("Error installing DirichletMultinomial")
    message(e$message)
    quit("no", -1)
  }
)

## Lastly, install all of the ArchR dependencies that arent installed by default:
library(ArchR)

#ArchR::installExtraPackages()
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

# We installed MACS2 using python and archr_requirements.txt
## It is also highly recommended that you install MACS2, which requires python,
## and have the macs2 executable in your PATH variable. This will allow ArchR to call peaks using MACS2.
##pip install MACS2
#
#devtools::install_github("GreenleafLab/chromVARmotifs")
tryCatch({
    devtools::install_github("GreenleafLab/chromVARmotifs")
},
    error = function(e) {
    message("Error installing GreenleafLab chromVARmotifs")
    message(e$message)
    quit("no", -1)
  }
)

#install.packages("magick")
tryCatch({
    install.packages("magick")
},
    error = function(e) {
    message("Error installing magick")
    message(e$message)
    quit("no", -1)
  }
)

# See section 1.5.1 Creating a Custom ArchGenome in
# https://www.archrproject.com/bookdown/getting-set-up.html
tryCatch({
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
},
    error = function(e) {
    message("Error installing NCBI GRCh38")
    message(e$message)
    quit("no", -1)
  }
)

tryCatch({
    BiocManager::install("org.Hs.eg.db")
},
    error = function(e) {
    message("Error installing org.Hs.eg.db")
    message(e$message)
    quit("no", -1)
  }
)

tryCatch({
    BiocManager::install("GenomicFeatures")
 },
    error = function(e) {
    message("Error installing GenomicFeatures")
    message(e$message)
    quit("no", -1)
  }
)   
