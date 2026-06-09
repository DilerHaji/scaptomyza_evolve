#!/bin/bash
# Post-deploy: install CRAN-only R packages + compile ngsLD
set -euo pipefail

echo "Installing LDheatmap from CRAN (not on conda)..."
Rscript -e 'install.packages("LDheatmap", repos="https://cloud.r-project.org")'

if [ -d "ngsLD-master" ] && [ ! -x "ngsLD-master/ngsLD" ]; then
    echo "Compiling ngsLD..."
    cd ngsLD-master
    make
    cd ..
    echo "ngsLD compiled: $(ls -la ngsLD-master/ngsLD 2>&1)"
fi

echo "Post-deploy complete."
