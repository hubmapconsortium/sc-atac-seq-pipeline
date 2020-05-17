FROM hubmap/sc-atac-seq

WORKDIR /opt/supplementary-data
RUN wget -O - https://giygas.compbio.cs.cmu.edu/sc-atac-seq/supplementary-data.tar.xz | tar -xJf -

WORKDIR /opt/supplementary-data/bwa-index
# compressing the index with xz reduces the archive size from 3.4GB to 2.9GB but
# requires 4-5x more time to decompress, so not really worth it
RUN wget -O - https://giygas.compbio.cs.cmu.edu/sc-atac-seq/grch38_index.tar.gz | tar -xzf -

WORKDIR /opt
