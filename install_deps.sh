#!/bin/bash
# One-time setup of the `crocotel` conda env + R deps. Idempotent: creates the
# env with R 4.4 if it does not exist, then installs the packages the crocotel
# R functions use. Run once on a fresh Hoffman2 account before the sweep.
set -eo pipefail
source ~/.bashrc

# 1. Create the env (R 4.4) if it is not already present.
if ! conda env list | awk '{print $1}' | grep -qx crocotel; then
  echo "Creating conda env 'crocotel' (R 4.4)..."
  conda create -n crocotel -c conda-forge -y r-base=4.4
fi

# 2. Install R package deps (bigsnpr pulls bigstatsr; data.table listed
#    explicitly since the sweep scripts use it directly).
conda install -n crocotel -c conda-forge -y \
  r-bigsnpr r-matrixeqtl r-glmnet r-data.table

echo "Done. Activate with: conda activate crocotel"
