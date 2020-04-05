FROM rocker/r-ver:3.6.3

LABEL maintainer="jshands@ucsc.edu"

# Use an env var so user is not requested to provide a
# time zone for tzdata
# https://stackoverflow.com/questions/44331836/apt-get-install-tzdata-noninteractive
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgsl-dev \
    liblzma-dev \
    libncurses-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    python-dev \
    python-pip \
    python3-dev \
    python3-pip \
    wget \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists

# Make sure a 'python' command is available so bedtools will install, as required in
# https://github.com/arq5x/bedtools2/blob/58e9973af1b3f5e3b26e5584aad7dc7b720f8765/Makefile#L192
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

COPY requirements_py2.txt /opt
COPY requirements_py3.txt /opt
RUN python2 -m pip install "Cython==0.29.15" --install-option="--no-cython-compile" \
 && python2 -m pip install -r /opt/requirements_py2.txt \
 # MACS2's `setup.py` file imports `numpy`, so need to install NumPy first.
 # See https://github.com/taoliu/MACS/issues/364
 && python3 -m pip install "numpy==1.18.2" \
 && python3 -m pip install -r /opt/requirements_py3.txt \
 && rm -rf /root/.cache/pip /opt/requirements_py2.txt /opt/requirements_py3.txt

COPY install_R_packages.R /opt
RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

WORKDIR /opt

ENV SAMTOOLS_VERSION=1.10
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && cd samtools-${SAMTOOLS_VERSION} \
 && ./configure \
 && make \
 && make install \
 && rm -rf /opt/samtools-${SAMTOOLS_VERSION} /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2

ENV BWA_VERSION=0.7.17
RUN wget -O bwa-${BWA_VERSION}.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download \
 && tar -xf bwa-${BWA_VERSION}.tar.bz2 \
 && cd bwa-${BWA_VERSION} \
 && make \
 && cp bwa /usr/local/bin \
 && rm -rf /opt/bwa-${BWA_VERSION} /opt/bwa-${BWA_VERSION}.tar.bz2

ENV BEDTOOLS_VERSION=2.29.1
# Install bedtools
RUN wget -O bedtools-${BEDTOOLS_VERSION}.tar.gz https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
 && tar -xf bedtools-${BEDTOOLS_VERSION}.tar.gz \
 && cd bedtools2 \
 && make \
 && make install \
 && rm -rf /opt/bedtools-${BEDTOOLS_VERSION}.tar.gz /opt/bedtools2

COPY create_genome_size_file.sh /tools/
COPY create_reference_genome_index.sh /tools/
COPY add_barcodes_to_reads.pl /tools/
COPY remove_blacklist.sh /tools/
COPY snapAnalysis_select_barcode.R /tools/
COPY snapAnalysis.R /tools/
COPY snapMotifAnalysis.R /tools/
COPY gather_sequence_files.py /tools/

ENV PATH=$PATH:/opt:/tools

CMD ["/bin/bash"]
