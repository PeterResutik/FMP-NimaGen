# Base image
FROM ubuntu:20.04

# Set non-interactive mode to prevent issues during package installation
ENV DEBIAN_FRONTEND=noninteractive

COPY scripts /workspace/scripts

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    unzip \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set up environment variables for Conda
ENV PATH="/opt/conda/bin:$PATH"

# Install Mamba (faster alternative to Conda)
RUN conda install -n base -c conda-forge mamba -y

# Create a Conda environment with specific versions of bioinformatics tools
RUN mamba create -n bioinfo -c bioconda -c conda-forge -y \
    python=3.9 \
    cutadapt=4.9 \
    samtools=1.21 \
    bwa=0.7.18

# Set default shell to activate Conda environment
SHELL ["/bin/bash", "-c"]

# Activate the bioinfo environment whenever the container is used
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bioinfo"]
CMD ["bash"]
