FROM rocker/r-ver:3.6.3

# Use an env var so user is not requested to provide a
# time zone for tzdata
# https://stackoverflow.com/questions/44331836/apt-get-install-tzdata-noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Install requirements
RUN apt update && apt install -y \
    build-essential \
    cython \
    gzip \
    htop \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libgsl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    pandoc \
    perl \
    python-pip \
    python3-pip \
    software-properties-common \
    texlive-latex-extra \
    unzip \
    wget \
    vim \
    zlib1g \
    zlib1g-dev


LABEL maintainer="jshands@ucsc.edu"


# Install libpng12 which is no longer available in Ubuntu package archive
# but is needed by other software that we will install
RUN wget https://launchpad.net/~ubuntu-security/+archive/ubuntu/ppa/+build/15108504/+files/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb 
RUN dpkg -i libpng12-0_1.2.54-1ubuntu1.1_amd64.deb

RUN apt install libreadline-dev

# Install R packages
RUN R -e "install.packages(c('Matrix', 'doSNOW', 'plot3D', 'optparse'))"

# Install SnapATAC
RUN R -e "install.packages(c('devtools'))"
RUN R -e "library(devtools); install_github('r3fang/SnapATAC')"

# Install R packages using BioConductor
# at https://bioconductor.org/install/#install-bioconductor-packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.10')"

# rtracklayer needed for cell selection work around described at
#      https://github.com/r3fang/SnapATAC/issues/139
# RFLPtools needed for the write.hclust call
# rmarkdown used when running an R markdown document
# dplyr needed for R markdown to create PDF files
RUN R -e "BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'rtracklayer', 'RFLPtools', 'rmarkdown', 'dplyr'))"

# For motif analysis
# Needed to install chromVAR
RUN R -e "BiocManager::install(c('motifmatchr'))"
RUN R -e "BiocManager::install(c('SummarizedExperiment'))"
# Needed to run chromVAR
RUN R -e "BiocManager::install(c('JASPAR2016'))"
RUN R -e "BiocManager::install(c('chromVAR'))"

RUN pip install html5lib

WORKDIR /opt/samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 -O samtools.tar.bz2 && \
    tar -xjvf samtools.tar.bz2 && \
    cd samtools-1.9 && \
    ./configure --prefix /opt/samtools/samtools-1.9 && \
    make && \
    make install
ENV PATH /opt/samtools/samtools-1.9/bin:$PATH

# Install BWA
WORKDIR /install
RUN wget -O "bwa-0.7.17.tar.bz2" "https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download" && \
    tar xvjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mkdir /tools/ && \
    cp bwa /tools/

# Install bedtools
WORKDIR /install
RUN wget -O "bedtools-2.29.1.tar.gz" https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
    tar -zxvf bedtools-2.29.1.tar.gz && \
    cd bedtools2 && \
    make && \
    cp bin/* /tools/


# Install SnapTools
# https://pypi.org/project/snaptools/
# Use pip instead of pip3 because snaptools requires python 2.7
RUN pip install 'snaptools==1.4.8'

# Python 3 required For MACS 2.2.6
# numpy needed for MACS2
RUN pip3 install numpy
RUN pip3 install 'MACS2==2.2.6'

COPY create_genome_size_file.sh /tools/
COPY create_reference_genome_index.sh /tools/
COPY add_barcodes_to_reads.pl /tools/
COPY remove_blacklist.sh /tools/
COPY snapAnalysis_select_barcode.R /tools/
COPY snapAnalysis.R /tools/
COPY snapMotifAnalysis.R /tools/
COPY gather_sequence_files.py /tools/

ENV PATH /tools/:$PATH



