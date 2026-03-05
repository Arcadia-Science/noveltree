# Use the official Python image as a base image
FROM python:3.9-slim

# Install procps for compatibility with nextflow
RUN apt-get update && \
    apt-get install -y \
        procps

# Install the required Python dependencies
# Install numpy first to avoid binary incompatibility issues
RUN pip install --no-cache-dir \
    numpy==1.24.3 \
    biopython==1.81 \
    pandas==2.0.0

# Copy the physicochemical properties script
COPY bin/zoogle/genefam_aa_summaries.py /usr/local/bin/genefam_aa_summaries.py
RUN chmod +x /usr/local/bin/genefam_aa_summaries.py

# Set the default command to launch an interactive shell
CMD ["/bin/bash", "-l"]