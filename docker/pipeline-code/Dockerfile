FROM hubmap/rocker-r-grch38:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libcairo2-dev \
    libcurl4-gnutls-dev \
    libncurses-dev \
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
#

RUN apt-get update \
 && apt-get -y install \
    libncurses-dev
#
WORKDIR /opt
#
ENV SAMTOOLS_VERSION=1.10
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
 && tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR samtools-${SAMTOOLS_VERSION}
RUN ./configure \
 && make \
 && make install \
 && rm -rf /opt/samtools-${SAMTOOLS_VERSION} /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2

WORKDIR /opt

ENV BWA_VERSION=0.7.17
RUN wget -O bwa-${BWA_VERSION}.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2/download \
 && tar -xf bwa-${BWA_VERSION}.tar.bz2
WORKDIR bwa-${BWA_VERSION}
RUN make \
 && cp bwa /usr/local/bin \
 && rm -rf /opt/bwa-${BWA_VERSION} /opt/bwa-${BWA_VERSION}.tar.bz2

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

# Get GTF for NCBI GRCh38 so we can make custom gene annotation for ArchR
# https://www.archrproject.com/bookdown/getting-set-up.html
# from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
#RUN wget -O GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz \
# && gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz

RUN wget -O GRCh37_latest_genomic.gff.gz https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz

COPY bin /opt

CMD ["/bin/bash"]
