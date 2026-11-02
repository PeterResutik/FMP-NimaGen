# The tools FMP-NimaGen runs, for -profile docker and -profile singularity (see README.md).
# The pipeline itself (main.nf, modules, resources) stays outside the image; Nextflow mounts it.
#
#   docker build --platform linux/amd64 -t peterresutik/nimagen-pipeline:<version> .
#
# Built for linux/amd64 only: rtn's SeqLib uses x86 SSE instructions. On Apple Silicon, Docker runs
# the image through emulation (slower, same results).
FROM mambaorg/micromamba:2.0.5@sha256:8d0d4cad074afe3ac93b58dce88bd2659153d3fa6de3a6968f501710e8d3c832

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        bzip2 procps ca-certificates git make autoconf automake g++ zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# bwa, samtools, FLASH, cutadapt, FDSTOOLS and Python with its packages: every package at the exact version
# of the lock files in containers/ (FMP-NimaGen.yml names only the main versions, for Conda on any platform)
COPY containers/FMP-NimaGen.linux-64.lock containers/FMP-NimaGen.pip.lock containers/rtn-build.linux-64.lock /tmp/
RUN micromamba create -y -n FMP-NimaGen -f /tmp/FMP-NimaGen.linux-64.lock \
    && /opt/conda/envs/FMP-NimaGen/bin/pip install --no-deps --no-cache-dir -r /tmp/FMP-NimaGen.pip.lock \
    && micromamba clean -a -y

# FDSTOOLS names alleles with STRNaming, which downloads the chrM reference from Ensembl the first
# time it needs it and keeps it in the home folder. A container may have neither internet access nor
# a home folder (Nextflow starts Singularity with --no-home), so chrM is downloaded here once and
# kept in STRNaming's own data folder.
RUN export HOME=/tmp/strnaming-home \
    && /opt/conda/envs/FMP-NimaGen/bin/strnaming refseq-cache chrM:1..16569 \
    && mv $HOME/.strnaming-cache/refseq-M-*.txt.gz \
       "$(/opt/conda/envs/FMP-NimaGen/bin/python -c 'import pathlib, strnaming; print(pathlib.Path(strnaming.__file__).parent / "data")')" \
    && rm -rf $HOME

# rtn with the three fixes (fork PeterResutik/RtN, branch fixes, with the upstream commits up to 4883c67),
# built for any x86-64 CPU; its source stays in /opt/RtN
RUN micromamba create -y -n rtn-build -f /tmp/rtn-build.linux-64.lock && micromamba clean -a -y \
    && git clone -b fixes https://github.com/PeterResutik/RtN.git /opt/RtN \
    && cd /opt/RtN && git checkout 1287f8d \
    && PREFIX=/opt/conda/envs/rtn-build MARCH=-march=x86-64-v2 ./build.sh \
    && cp rtn /usr/local/bin/rtn

ENV PATH=/opt/conda/envs/FMP-NimaGen/bin:$PATH
