FROM ubuntu

RUN	apt-get -y update
RUN mkdir -p ~/src

# Install 'bwa'

RUN apt-get install -y \
  git \
  build-essential \
  zlib1g-dev

RUN	git clone https://github.com/lh3/bwa.git ~/src/bwa
RUN make -C ~/src/bwa


# Install 'samtools'

RUN apt-get install -y \
  libncurses5-dev \
  libncursesw5-dev \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  autoconf

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
