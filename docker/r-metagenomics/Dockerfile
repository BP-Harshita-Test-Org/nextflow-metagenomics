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
), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/'); \
    BiocManager::install(c( \
    'phyloseq', \
    'ANCOMBC', \
    'Maaslin2', \
    'ComplexHeatmap' \
), ask=FALSE, update=FALSE)"

WORKDIR /data
