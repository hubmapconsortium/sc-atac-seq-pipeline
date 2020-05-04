FROM hubmap/sc-atac-seq

WORKDIR /opt
RUN curl https://giygas.compbio.cs.cmu.edu/sc-atac-seq/prereqs.tar | tar -xf -

WORKDIR /opt/bwa-index
RUN curl https://giygas.compbio.cs.cmu.edu/sc-atac-seq/grch38_index.tar.gz | tar -xzf -
