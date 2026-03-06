FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        git \
        wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install PREQUAL
RUN git clone https://github.com/simonwhelan/prequal.git /opt/prequal && \
    cd /opt/prequal && \
    make && \
    cp prequal /usr/local/bin/ && \
    rm -rf /opt/prequal

CMD ["/bin/bash", "-l"]
