#!/bin/bash

set -e

echo "Installing system dependencies..."

apt-get update

apt-get install -y \
    build-essential \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libuv1-dev


echo "Installing R packages..."

R -e "install.packages(c(
    'tidyverse',
    'languageserver',
    'here',
    'renv',
    'BiocManager'
))"


echo "Setup complete."
