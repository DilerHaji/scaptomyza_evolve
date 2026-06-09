#!/bin/bash
# Post-deploy script for the haplovalidate conda environment.
# Snakemake runs this automatically after creating the env from haplovalidate.yaml.
# Installs the three non-CRAN R packages from local source directories.

set -euo pipefail

# Navigate to the project root (Snakemake runs post-deploy from the working dir)
echo "Installing local R packages into conda env..."

# haploReconstruct first (dependency of haplovalidate)
if [ -d "haploReconstruct-master" ]; then
    Rscript -e 'install.packages("haploReconstruct-master", repos=NULL, type="source")'
    echo "  haploReconstruct installed"
fi

# ACER
if [ -d "ACER-master" ]; then
    Rscript -e 'install.packages("ACER-master", repos=NULL, type="source")'
    echo "  ACER installed"
fi

# haplovalidate
if [ -d "haplovalidate-master" ]; then
    Rscript -e 'install.packages("haplovalidate-master", repos=NULL, type="source")'
    echo "  haplovalidate installed"
fi

echo "Post-deploy complete."
