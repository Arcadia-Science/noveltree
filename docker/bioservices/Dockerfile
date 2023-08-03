# Use the official Python image as a base image
FROM python:3.9-slim

# Install procps for compatibility with nextflow 
RUN apt-get update && \
    apt-get install -y \
        procps
        
# Install the required Python dependencies
RUN pip install --no-cache-dir \
    argparse=1.1 \
    bioservices==1.10.0 \
    pandas==2.0.0

# Set the default command to launch an interactive shell
CMD ["/bin/bash", "-l"]
