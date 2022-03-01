FROM rocker/tidyverse:4.1.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libcairo2-dev \
    libcurl4-gnutls-dev \
    libncurses-dev \
    libglpk-dev \
    libgsl-dev \
    libmagick++-dev \
    libxml2-dev \
    libssl-dev \
    python3-dev \
    python3-pip \
 && rm -rf /var/lib/apt/lists

# Make sure a 'python' command is available so bedtools will install, as required in
# https://github.com/arq5x/bedtools2/blob/58e9973af1b3f5e3b26e5584aad7dc7b720f8765/Makefile#L192
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

COPY requirements.txt /opt
RUN python3 -m pip install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip /opt/requirements.txt

RUN apt-get update \
 && apt-get -y install \
    libncurses-dev

WORKDIR /opt
ENV SAMTOOLS_VERSION=1.13
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR samtools-${SAMTOOLS_VERSION}
RUN ./configure \
 && make \
 && make install \
 && rm -rf /opt/samtools-${SAMTOOLS_VERSION} /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2

WORKDIR /opt
ENV BEDTOOLS_VERSION=2.29.1
# Install bedtools
RUN wget -O bedtools-${BEDTOOLS_VERSION}.tar.gz https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz \
 && tar -xf bedtools-${BEDTOOLS_VERSION}.tar.gz
WORKDIR bedtools2
RUN make \
 && make install \
 && rm -rf /opt/bedtools-${BEDTOOLS_VERSION}.tar.gz /opt/bedtools2

WORKDIR /opt
ENV HTSLIB_VERSION=1.13
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
 && tar -xf htslib-${HTSLIB_VERSION}.tar.bz2
WORKDIR htslib-${HTSLIB_VERSION}
RUN ./configure \
 && make \
 && make install \
 && rm -rf /opt/htslib-${HTSLIB_VERSION} /opt/htslib-${HTSLIB_VERSION}.tar.bz2

WORKDIR /opt

RUN wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gff3.gz

COPY bin /opt

CMD ["/bin/bash"]
