# Base image
FROM ubuntu:20.04

# Set non-interactive mode
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    bzip2 \
    ca-certificates \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set up Conda path
ENV PATH="/opt/conda/bin:$PATH"

# Install Mamba
RUN conda install -n base -c conda-forge mamba -y

# Create Conda environment
RUN mamba create -n bioinfo -c bioconda -c conda-forge -y \
    python=3.9 \
    cutadapt=4.6 \
    samtools=1.21 \
    bwa=0.7.18

# Copy scripts into the container
COPY scripts /workspace/scripts

# Set default shell
SHELL ["/bin/bash", "-c"]

# Activate environment when container runs
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bioinfo"]
CMD ["bash"]

