FROM nciccbr/ccbr_ubuntu_22.04:v4
#FROM continuumio/conda-ci-linux-64-python3.9  

ENV TAG="v0.1"
WORKDIR /opt2

# Set environment variables
ENV JULIA_VERSION=1.6.7
ENV JULIA_DIR=/opt/julia

# Install miniconda
ENV CONDA_DIR /opt/conda

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
     /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH


RUN conda config  --append channels conda-forge && \
    conda config  --append channels bioconda && \
    conda config  --append channels conda-forge && \
    conda config  --set restore_free_channel true


#RUN conda install conda-forge::julia-1.7.1

#RUN echo | curl -fsSL https://install.julialang.org | sh
#RUN juliaup add 1.6

# Download and extract Julia
RUN wget -q https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.7-linux-x86_64.tar.gz && \
    tar -xzf julia-$JULIA_VERSION-linux-x86_64.tar.gz && \
    mv julia-$JULIA_VERSION /opt/julia && \
    ln -s /opt/julia/bin/julia /usr/local/bin/julia && \
    rm julia-$JULIA_VERSION-linux-x86_64.tar.gz

RUN julia --version

# Install Whippet
RUN git clone https://github.com/timbitz/Whippet.jl.git
RUN cd Whippet.jl
RUN julia +1.6 --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'



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
