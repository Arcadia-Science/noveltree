FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir orthoxml-tools==1.3.0 lxml
CMD ["/bin/bash", "-l"]
