FROM rocker/r-ver:3.6.3

#FROM rocker/r-ver:4.0.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    curl \
    xz-utils \
 && rm -rf /var/lib/apt/lists

WORKDIR /opt/supplementary-data
RUN curl https://giygas.compbio.cs.cmu.edu/sc-atac-seq/supplementary-data.tar.xz | tar -xJf -

WORKDIR /opt/supplementary-data/bwa-index
# compressing the index with xz reduces the archive size from 3.4GB to 2.9GB but
# requires 4-5x more time to decompress, so not really worth it
RUN curl https://giygas.compbio.cs.cmu.edu/sc-atac-seq/grch38_index.tar.gz | tar -xzf -

WORKDIR /opt
