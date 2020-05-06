FROM hubmap/sc-atac-seq

WORKDIR /opt/supplementary-data
RUN wget -O - https://giygas.compbio.cs.cmu.edu/sc-atac-seq/supplementary-data.tar.xz | tar -xJf -

WORKDIR /opt/supplementary-data/bwa-index
RUN wget -O - https://giygas.compbio.cs.cmu.edu/sc-atac-seq/grch38_index.tar.gz | tar -xzf -

WORKDIR /opt
