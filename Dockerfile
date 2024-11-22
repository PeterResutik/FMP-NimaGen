# Start with a basic Ubuntu image
FROM ubuntu:20.04

# Set environment variable to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install bwa and samtools
RUN apt-get update && \
    apt-get install -y bwa samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to bash
CMD ["bash"]
