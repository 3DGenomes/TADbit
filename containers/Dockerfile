# Use Debian Bullseye as base image
FROM debian:bullseye-slim

# Set environmental variables
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/usr/local/anaconda/envs/env/bin:$PATH

# Install necessary packages for setting up the environment
RUN apt-get update -q && \
    apt-get install -q -y --no-install-recommends \
        dirmngr \
        apt-transport-https \
        ca-certificates \
        software-properties-common \
        gnupg2 \
        wget \
        bzip2 \
        unzip \
        build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget --quiet --no-check-certificate "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" \
      -O /root/miniconda.sh && \
    /bin/bash /root/miniconda.sh -b -p /usr/local/anaconda && \
    rm /root/miniconda.sh

# Initialize Conda and update it
RUN /usr/local/anaconda/bin/conda init bash && \
    /usr/local/anaconda/bin/conda config --set always_yes yes --set changeps1 no && \
    /usr/local/anaconda/bin/conda update -q conda

# Add Conda channels
RUN /usr/local/anaconda/bin/conda config --add channels conda-forge && \
    /usr/local/anaconda/bin/conda config --add channels bioconda && \
    /usr/local/anaconda/bin/conda config --add channels salilab

# Create and activate the Conda environment
# Install necessary packages within the environment
RUN /usr/local/anaconda/bin/conda create -n env python=3.9 && \
    /usr/local/anaconda/bin/conda install -n env -y -q \
        imp \
        scipy \
        numpy \
        matplotlib \
        mcl \
        samtools \
        pysam \
        ucsc-bigWigToBedGraph \
        ucsc-wigtobigwig \
        gem3-mapper \
        hisat2

# Clean up Conda to reduce image size
RUN /usr/local/anaconda/bin/conda clean -y --all

# Final touches and setting up the environment for when the container runs
WORKDIR /data
ENTRYPOINT ["/usr/local/anaconda/envs/env/bin/conda", "run", "-n", "env"]
CMD ["/bin/bash"]

