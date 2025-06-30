FROM condaforge/miniforge3:latest

ENV TAG="v0.1"
WORKDIR /opt2
RUN cd /opt2

RUN mamba install conda-forge::r-base=4.2.1 conda-forge::r-rcpp=1.0.10 conda-forge::r-rstan=2.21.7 conda-forge::r-stanheaders=2.21.0-7 conda-forge::r-devtools conda-forge::r-hmisc conda-forge::r-reshape2 conda-forge::r-foreach conda-forge::r-dplyr conda-forge::r-domc conda-forge::r-shinyjs r::r-dt conda-forge::intervals conda-forge::r-shinycssloaders conda-forge::r-bh

RUN conda install conda-forge::r-r.utils

RUN Rscript -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages("optparse", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("Biobase")'
RUN Rscript -e 'install.packages("TailRank", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'BiocManager::install("DirichletMultinomial")'

RUN conda install conda-forge::r-devtools
RUN Rscript -e 'devtools::install_github("stan-dev/rstantools")'
RUN Rscript -e 'install.packages("intervals", repos = "http://cran.us.r-project.org")'
RUN Rscript -e 'devtools::install_github("davidaknowles/leafcutter/leafcutter", ref = "psi_2019")'

RUN git clone https://github.com/davidaknowles/leafcutter


# cleanup etc
COPY Dockerfile.v0.1 /opt2/Dockerfile.v0.1
RUN chmod a+rX /opt2/Dockerfile.v0.1
WORKDIR /data2
RUN apt-get clean

# details
LABEL maintainer="guibletwm@nih.gov"
LABEL github_handle="wilfriedguiblet"
ARG BUILD_DATE
LABEL build_date=$BUILD_DATE