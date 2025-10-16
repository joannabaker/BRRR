# Dev container with R + system deps for Shiny and your packages
FROM mcr.microsoft.com/devcontainers/base:ubuntu

# System libs for tidyverse, DT, curl/ssl, XML, graphics, git2, etc.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
      r-base r-base-dev \
      libcurl4-openssl-dev libssl-dev libxml2-dev \
      libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
      libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
      libcairo2-dev libxt6 libgit2-dev pandoc && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Pre-create an R library path inside the container (optional)
ENV R_LIBS_USER=/usr/local/lib/R/site-library
