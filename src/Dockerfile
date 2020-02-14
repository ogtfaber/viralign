FROM ubuntu


# Setup

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

RUN mkdir -p ~/bin


# Install 'bwa'

RUN	git clone https://github.com/lh3/bwa.git ~/bin/bwa
RUN make -C ~/bin/bwa


# Install 'samtools'

RUN	git clone https://github.com/samtools/htslib ~/bin/htslib

WORKDIR /root/bin/htslib
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install

RUN	git clone https://github.com/samtools/samtools ~/bin/samtools
WORKDIR /root/bin/samtools
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install


# Download and index virus reference file

# RUN mkdir -p ~/data
# WORKDIR /root/data
# RUN wget https://rvdb.dbi.udel.edu/download/U-RVDBv17.0.fasta.gz
# RUN ~/bin/bwa/bwa index ~/data/U-RVDBv17.0.fasta.gz
# RUN mkdir -p ~/ref
# COPY ref /root


# Run the pipeline

COPY run.sh /root/run.sh
RUN chmod +x ~/run.sh

CMD ["/root/run.sh"]

# docker run \
#   -v /Users/onnofaber/Projects/rm/viralign/sample:/data \
#   -v /Users/onnofaber/Projects/rm/viralign/ref:/ref \
#   my-image
