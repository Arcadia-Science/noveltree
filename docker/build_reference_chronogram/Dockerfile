# Use the official Python image as a base image
FROM python:3.9-slim

# Install procps for compatibility with nextflow, and gcc for biopython build
RUN apt-get update && \
    apt-get install -y \
        procps \
        gcc \
        g++ && \
    rm -rf /var/lib/apt/lists/*

# Install the required Python dependencies
RUN pip install --upgrade pip && \
    pip install --no-cache-dir \
    numpy==1.24.3 \
    scipy==1.10.1 \
    biopython==1.83

# Copy the build_reference_chronogram script
COPY bin/build_reference_chronogram.py /usr/local/bin/build_reference_chronogram.py
RUN chmod +x /usr/local/bin/build_reference_chronogram.py

# Set the default command to launch an interactive shell
CMD ["/bin/bash", "-l"]
