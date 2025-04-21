FROM nciccbr/ccbr_ubuntu_22.04:v4
#FROM continuumio/conda-ci-linux-64-python3.9  

ENV TAG="v0.1"
WORKDIR /opt2


# Install miniconda
ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH


RUN conda config  --append channels conda-forge && \
    conda config  --append channels bioconda && \
    conda config  --append channels conda-forge && \
    conda config  --set restore_free_channel true


RUN conda install conda-forge::julia-1.7.1

# Install Whippet
RUN git clone https://github.com/timbitz/Whippet.jl.git && \
    cd Whippet.jl && \
    julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'



# cleanup etc
COPY Dockerfile /opt2/Dockerfile
RUN chmod a+rX /opt2/Dockerfile
WORKDIR /data2
RUN apt-get clean

# details
LABEL maintainer="guibletwm@nih.gov"
LABEL github_handle="wilfriedguiblet"
ARG BUILD_DATE
LABEL build_date=$BUILD_DATE
