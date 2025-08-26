# Use official slim Python 3.12 image
FROM python:3.12-slim

# Install system build dependencies needed for scientific Python stack
RUN apt-get update && apt-get install -y \
      build-essential \
      python3-dev \
      libxml2-dev \
      libxslt-dev \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libglpk40 \
      libffi-dev \
      libssl-dev \
      libhdf5-dev hdf5-tools \
      pkg-config \
      git \
      procps \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install Python requirements
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade pip setuptools wheel \
    && pip install -r /tmp/requirements.txt

# Working directory inside the container
WORKDIR /workspace

# Default command (optional)
CMD ["bash"]