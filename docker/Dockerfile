FROM ubuntu

RUN apt-get update && apt-get install -y \
    git \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    autoconf \
    wget \
  && rm -rf /var/lib/apt/lists/*

RUN mkdir -p ~/src


# Install 'bwa'
RUN	git clone https://github.com/lh3/bwa.git ~/src/bwa
RUN make -C ~/src/bwa

# Install 'samtools'
RUN	git clone https://github.com/samtools/htslib ~/src/htslib

WORKDIR /root/src/htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install

RUN	git clone https://github.com/samtools/samtools ~/src/samtools
WORKDIR /root/src/samtools
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install


# Download virus reference file

RUN mkdir -p ~/data
WORKDIR /root/data
RUN wget https://rvdb.dbi.udel.edu/download/U-RVDBv17.0.fasta.gz

# RUN ~/src/bwa/bwa index ~/data/U-RVDBv17.0.fasta.gz

# CMD ["ls -al"]