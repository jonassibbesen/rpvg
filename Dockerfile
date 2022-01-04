FROM mirror.gcr.io/library/ubuntu:20.04

MAINTAINER jsibbese@ucsc.edu

WORKDIR /home

### Install essential tools (including protobuf)

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends wget git autoconf automake make cmake gcc g++ pkg-config protobuf-compiler libprotoc-dev libprotobuf-dev libjansson-dev && \
    rm -rf /var/lib/apt/lists/*

### Install htslib dependencies

RUN apt-get update && \
	DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev && \
    rm -rf /var/lib/apt/lists/*

### Compile rpvg 

RUN GIT_SSL_NO_VERIFY=true git clone --recursive https://github.com/jonassibbesen/rpvg.git && \
	cd rpvg && \
	mkdir build && \
	cd build && \
	cmake .. && \
	make && \
	cd ../../ && \
	mv rpvg/bin/rpvg /usr/bin/ && \
	unlink rpvg/deps/libvgio/vg
