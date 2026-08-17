#!/bin/bash
# One-time install of R deps into the crocotel conda env.

source ~/.bashrc
conda install -n crocotel -c conda-forge -y \
  r-bigsnpr r-matrixeqtl r-glmnet
