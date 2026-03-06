FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        wget \
        zlib1g-dev \
        python3 \
        python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install SeqKit v2.6.1
RUN wget -q https://github.com/shenwei356/seqkit/releases/download/v2.6.1/seqkit_linux_amd64.tar.gz && \
    tar -xzf seqkit_linux_amd64.tar.gz && \
    mv seqkit /usr/local/bin/ && \
    rm seqkit_linux_amd64.tar.gz

# Install CD-HIT v4.8.1
RUN wget -q https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz && \
    tar -xzf cd-hit-v4.8.1-2019-0228.tar.gz && \
    cd cd-hit-v4.8.1-2019-0228 && \
    make && \
    cp cd-hit cd-hit-est /usr/local/bin/ && \
    cd / && \
    rm -rf cd-hit-v4.8.1-2019-0228*

# Install TransDecoder v5.7.1
RUN wget -q https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz && \
    tar -xzf TransDecoder-v5.7.1.tar.gz && \
    mv TransDecoder-TransDecoder-v5.7.1 /opt/TransDecoder && \
    rm TransDecoder-v5.7.1.tar.gz && \
    ln -s /opt/TransDecoder/TransDecoder.LongOrfs /usr/local/bin/ && \
    ln -s /opt/TransDecoder/TransDecoder.Predict /usr/local/bin/

CMD ["/bin/bash", "-l"]
