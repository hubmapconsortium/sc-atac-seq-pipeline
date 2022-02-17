FROM pypy:3

WORKDIR /opt/supplementary-data/hisat2-index
RUN curl -s https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz | tar -xzf -

WORKDIR /opt

RUN curl -o hisat2.zip https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
 && unzip hisat2.zip \
 && rm hisat2.zip

ENV PATH="/opt/hisat2-2.2.1:${PATH}"

ENV SAMTOOLS_VERSION=1.10
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && cd samtools-${SAMTOOLS_VERSION} \
 && ./configure \
 && make \
 && make install \
 && rm -rf /opt/samtools-${SAMTOOLS_VERSION} /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2

COPY bin /opt

CMD ["/bin/bash"]
