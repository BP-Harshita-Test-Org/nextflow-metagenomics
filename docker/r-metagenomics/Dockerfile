FROM rocker/r-ver:4.3.2

LABEL maintainer="ASD Biomarker Study Pipeline"
LABEL description="R environment with metagenomics packages for diversity, MaAsLin2, and visualization"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Use Posit Package Manager for pre-compiled R binaries
ENV CRAN_REPO="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

# Install CRAN packages
RUN R -e "install.packages(c( \
    'optparse', \
    'dplyr', \
    'tidyr', \
    'readr', \
    'ggplot2', \
    'pheatmap', \
    'reshape2', \
    'vegan' \
), repos=Sys.getenv('CRAN_REPO'))"

# Install Bioconductor packages (split for better Docker caching)
RUN R -e "install.packages('BiocManager', repos=Sys.getenv('CRAN_REPO'))"
RUN R -e "BiocManager::install('phyloseq', ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install('ANCOMBC', ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install('Maaslin2', ask=FALSE, update=FALSE)"
RUN R -e "BiocManager::install('ComplexHeatmap', ask=FALSE, update=FALSE)"

WORKDIR /data
