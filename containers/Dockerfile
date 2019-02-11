FROM debian:8

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 --no-install-recommends && \
    rm -rf /var/lib/apt/lists/*

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet --no-check-certificate https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh \
        -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda config --add channels salilab && conda config --add channels bioconda && \
    conda install -y -q imp scipy matplotlib jupyter mcl samtools sra-tools pysam && \
    conda clean -y --all  && rm -rf /opt/conda/pkgs/*

RUN wget --quiet --no-check-certificate https://newcontinuum.dl.sourceforge.net/project/gemlibrary/gem-library/Binary%20pre-release%202/GEM-binaries-Linux-x86_64-core_i3-20121106-022124.tbz2 \
        -O GEM.tbz2 && \
    tar xvf GEM.tbz2 && cd GEM-*/ && \
    mv * /usr/local/bin/ && cd .. && rm -rf GEM*

RUN apt-get update --fix-missing && \
    apt-get install -y unzip build-essential --no-install-recommends && \
    wget --quiet --no-check-certificate https://github.com/fransua/TADbit/archive/dev.zip && unzip dev.zip && \
    cd TADbit-dev && python setup.py install && cd .. && rm -rf TADbit-dev dev.zip && \
    apt-get remove -y --purge unzip build-essential && \
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD [ "/bin/bash" ]
