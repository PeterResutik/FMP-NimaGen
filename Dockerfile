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
    build-essential \
    cmake \
    bc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set up Conda path
ENV PATH="/opt/conda/bin:$PATH"

# Install Mamba
RUN conda install -n base -c conda-forge mamba -y

# Create Conda environment with additional packages
RUN mamba create -n bioinfo -c bioconda -c conda-forge -y \
    python=3.9 \
    cutadapt=4.6 \
    samtools=1.21 \
    bwa=0.7.18 \
    flash=1.2.11 \
    csvtk \
    fastqc \
    pandas \
    biopython \
    matplotlib \
    openpyxl \
    gatk4 \
    bcftools

# Clone RtN repository and build it
RUN git clone https://github.com/Ahhgust/RtN.git /workspace/RtN 

# Copy the precompiled RtN binary from Nix_binary to a system-wide location
RUN cp /workspace/RtN/Nix_binary/rtn /usr/local/bin/ && chmod +x /usr/local/bin/rtn


# Set environment variable for RtN
ENV RTN_PATH="/workspace/RtN/build"

# Activate the environment and install fdstools via pip
RUN /opt/conda/envs/bioinfo/bin/pip install fdstools

# Copy scripts into the container
COPY resources/scripts /workspace/resources/scripts
COPY resources/primers /workspace/resources/primers
COPY resources/rtn_files /workspace/resources/rtn_files
COPY resources/amplicon_bed /workspace/resources/amplicon_bed
COPY resources/rCRS /workspace/resources/rCRS
COPY resources/fdstools /workspace/resources/fdstools

# Set default shell
SHELL ["/bin/bash", "-c"]

# Activate environment when container runs
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "bioinfo"]
CMD ["bash"]

