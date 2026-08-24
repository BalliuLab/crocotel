#!/bin/bash
# R CMD check --as-cran of the crocotel package on the cluster (release gate).
#$ -N crocotel_check
#$ -cwd
#$ -j y
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 4
source ~/.bashrc
conda activate crocotel
set -eo pipefail
export OMP_NUM_THREADS=${NSLOTS:-1}
export OPENBLAS_NUM_THREADS=${NSLOTS:-1}
: "${CROCOTEL_REPO:=/u/home/b/bballiu/balliulab/bballiu/crocotel}"
WORK=/u/scratch/${USER:0:1}/${USER}/pkg_check
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
R CMD build "$CROCOTEL_REPO/r_package"
R CMD check --as-cran --no-manual crocotel_*.tar.gz   # cluster TeX lacks inconsolata; manual is checked locally
tail -5 crocotel.Rcheck/00check.log
