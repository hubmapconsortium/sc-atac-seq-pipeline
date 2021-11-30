FROM hubmap/rocker-r-grch38:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    fastqc \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    libgsl-dev \
    liblzma-dev \
    libncurses-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    pkg-config \
    python3-dev \
    python3-pip \
    tabix \
    tk \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists

# Make sure a 'python' command is available so bedtools will install, as required in
# https://github.com/arq5x/bedtools2/blob/58e9973af1b3f5e3b26e5584aad7dc7b720f8765/Makefile#L192
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

COPY requirements.txt /opt
RUN python3 -m pip install "Cython==0.29.17" --install-option="--no-cython-compile" \
 && python3 -m pip install numpy \
 && python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip /opt/requirements.txt

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
RUN wget https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2 \
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

COPY bin /opt

CMD ["/bin/bash"]
