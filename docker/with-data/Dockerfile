FROM hubmap/sc-atac-seq

WORKDIR /opt
RUN wget -O - https://giygas.compbio.cs.cmu.edu/sc-atac-seq/supplementary-data.tar.xz | tar -xJf -

WORKDIR /opt/bwa-index
RUN wget -O - https://giygas.compbio.cs.cmu.edu/sc-atac-seq/grch38_index.tar.gz | tar -xzf -
