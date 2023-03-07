FROM mirror.gcr.io/library/ubuntu:20.04

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install essential tools (including protobuf and jansson)

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends wget git autoconf automake make cmake gcc g++ pkg-config protobuf-compiler libprotoc-dev libprotobuf-dev libjansson-dev && \
    rm -rf /var/lib/apt/lists/*

### Install htslib

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
	tar -xvjf htslib-1.16.tar.bz2 && \
	cd htslib-1.16 && \
	./configure && \
	make install && \
 	cd .. && \
 	rm -r *

### Get rpvg 

RUN GIT_SSL_NO_VERIFY=true git clone --recursive https://github.com/jonassibbesen/rpvg.git

### Set target CPU architecture

RUN cd rpvg && \
	sed -i 's/march=native/march=nehalem/g' CMakeLists.txt deps/sdsl-lite/CMakeLists.txt

### Compile rpvg 

RUN cd rpvg && \
	mkdir build && \
	cd build && \
	cmake .. && \
	make && \
	cd ../../ && \
	mv rpvg/bin/rpvg /usr/bin/ && \
	unlink rpvg/deps/libvgio/vg
