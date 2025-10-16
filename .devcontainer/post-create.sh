#!/usr/bin/env bash
set -e

# Install required R packages once container is built.
Rscript -e 'install.packages(c(
  "shiny","ape","ratematrix","ggplot2","dplyr","tidyr","purrr","DT",
  "colourpicker","readr","remotes"
), repos="https://cloud.r-project.org")'

# If BayesTraitR is GitHub-only, install like this (comment out if on CRAN):
Rscript -e 'if (!requireNamespace("BayesTraitR", quietly=TRUE)) remotes::install_github("mrborges23/BayesTraitR")'
chmod +x .devcontainer/post-create.sh
