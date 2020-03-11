FROM ubuntu:18.04

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
    libcurl3 \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libgsl-dev \
    libssl-dev \
    libxml2-dev \
    pandoc \
    perl \
    python3-pip \
    software-properties-common \
    texlive-latex-extra \
    unzip \
    wget \
    vim \
    zlib1g \
    zlib1g-dev


LABEL maintainer="jshands@ucsc.edu"

# Need to have both libcurl3 and libcurl4
# Move libcurl3 to a different location so we can install libcurl4
# https://dev.to/jake/using-libcurl3-and-libcurl4-on-ubuntu-1804-bionic-184g
RUN cp /usr/lib/x86_64-linux-gnu/libcurl.so.4.5.0 /usr/lib/libcurl.so.3
# Now install libcurl4
RUN apt install libcurl4 libcurl4-openssl-dev -y

# Make python3 the default so latest MACS version will install
# https://stackoverflow.com/questions/41986507/unable-to-set-default-python-version-to-python3-in-ubuntu/41986843 
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

# Install R 3.5 so we can use newer version of BioConductor
# https://bioconductor.org/install/ 
# https://www.r-bloggers.com/installation-of-r-3-5-on-ubuntu-18-04-lts-and-tips-for-spatial-packages/
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

# Install libpng12 which is no longer available in Ubuntu package archive
# but is needed by other software that we will install
RUN wget https://launchpad.net/~ubuntu-security/+archive/ubuntu/ppa/+build/15108504/+files/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb 
RUN dpkg -i libpng12-0_1.2.54-1ubuntu1.1_amd64.deb

# Install R 
# https://www.r-bloggers.com/installation-of-r-3-5-on-ubuntu-18-04-lts-and-tips-for-spatial-packages/
RUN apt install libreadline-dev
RUN apt install -y r-base-core r-base-dev r-base r-recommended

WORKDIR install

RUN pip3 install Cython --install-option="--no-cython-compile"

# Install SnapTools
# https://pypi.org/project/snaptools/
RUN pip3 install 'snaptools==1.4.8'

RUN pip3 install 'MACS2==2.2.6'

# Install R packages
RUN R -e "install.packages(c('devtools', 'Matrix', 'doSNOW', 'plot3D', 'optparse'))"
# Install SnapATAC
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

RUN pip3 install html5lib

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

COPY create_genome_size_file.sh /tools/
COPY create_reference_genome_index.sh /tools/
COPY add_barcodes_to_reads.pl /tools/
COPY remove_blacklist.sh /tools/
COPY snapAnalysis_select_barcode.R /tools/
COPY snapAnalysis.R /tools/

ENV PATH /tools/:$PATH
