# Use the official Python image as a base image
FROM python:3.9-slim

# Install procps for compatibility with nextflow 
RUN apt-get update && \
    apt-get install -y \
        procps
        
# Install the required Python dependencies
RUN pip install --no-cache-dir \
    argparse \
    bioservices==1.10.0 \
    pandas

# Set the default command to launch an interactive shell
CMD ["/bin/bash", "-l"]