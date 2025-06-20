FROM nciccbr/ccbr_ubuntu_22.04:v4

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


RUN Rscript -e 'install.packages("stan", repos = "http://cran.us.r-project.org") '&& \
    Rscript -e 'devtools::install_github("stan-dev/rstantools")' && \
    Rscript -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter")'

RUN conda install bioconductor-biobase bioconductor-dirichletmultinomial r-rcpp=1.0.10 r-rstan=2.21.7 r-stanheaders=2.21.0-7

 





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
