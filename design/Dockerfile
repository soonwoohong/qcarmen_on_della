# Latch preamble (do not change)
FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:ace9-main

ENV \
    WORK_DIR=/adapt \
    MEMO_DIR=/memo \
    MAFFT_PATH=/usr/bin/mafft

# Install base utilities
RUN apt-get update && \
    apt-get install -y wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# # The series of commands below installs 'Miniconda', which provides conda
# ENV CONDA_DIR /opt/conda
# ENV PATH=$CONDA_DIR/bin:$PATH

# Install base utilitiesRUN apt-get update && \apt-get install -y wget && \apt-get clean && \rm -rf /var/lib/apt/lists/*
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && mkdir /root/.conda && bash Miniconda3-latest-Linux-x86_64.sh -b && rm -f Miniconda3-latest-Linux-x86_64.sh
RUN /root/miniconda3/bin/conda create -y -n enrich2 -c defaults python=3.8
ENV PATH="/root/miniconda3/envs/enrich2/bin:${PATH}"

RUN /root/miniconda3/bin/conda install -n enrich2 numpy=1.18.5
RUN /root/miniconda3/bin/conda install -n enrich2 -c bioconda adapt=1.4.1
RUN /root/miniconda3/bin/conda install -n enrich2 -c conda-forge biopython=1.79
RUN /root/miniconda3/bin/conda install -n enrich2 -c bioconda primer3-py

# Install MAFFT
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.487-1_amd64.deb \
    && dpkg -i mafft_7.487-1_amd64.deb \
    && rm -rf mafft_7.487-1_amd64.deb

# STOP HERE:
# The following lines are needed to ensure your build environment works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root