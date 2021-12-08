FROM alpine:3.13.1

# R + DESeq2 + FactoMineR + factoextra + pheatmap + RColorBrewer + dplyr +xlsx
# https://cran.r-project.org/
# v4.0.3 / 3.3.2 / 2.3 / 1.0.7 / 1.0.12 / 1.1-2 / 1.0.2

ENV R_PKGS R=~4.0.3 \
  R-mathlib=~4.0.3 \
  R-dev=~4.0.3 \
  R-doc=~4.0.3

ENV R_DEPS g++ \
  libxml2-dev \
  make \
  cmake \
  linux-headers \
  cairo-dev \
  libxmu-dev \
  pango-dev \
  perl \
  tiff-dev \
  icu-dev \
  libjpeg-turbo \
  pcre-dev \
  readline-dev \
  libexecinfo-dev \
  file \
  ttf-linux-libertine \
  git \
  openssl \
  autoconf \
  automake \
  libuv \
  http-parser \
  tzdata \
  libgit2-dev \
  bash

RUN  apt-get update \
    && apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl git libopenmpi-dev openmpi-bin openmpi-doc libxml2-dev \
    && R -e "if (!requireNamespace('BiocManager', quietly = TRUE)){install.packages('BiocManager')}" \
    && R -e "BiocManager::install('DESeq2')" \
    && R -e "url='https://cran.r-project.org/src/contrib/pheatmap_1.0.12.tar.gz';install.packages(url, repos=NULL, type='source')" \ 
    && R -e "url='https://cran.r-project.org/src/contrib/RColorBrewer_1.1-2.tar.gz';install.packages(url, repos=NULL, type='source')" \
    && R -e "BiocManager::install('factoextra')" \
    && R -e "BiocManager::install('FactoMineR')" \
    && R -e "BiocManager::install('xlsx')"
    
ENTRYPOINT ["/usr/src/.docker_modules/r-base/4.0.3"]   
