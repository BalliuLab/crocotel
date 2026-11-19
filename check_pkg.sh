#!/bin/bash
# R CMD check --as-cran of the crocotel package on Hoffman2 (release gate).
#
# Submit from the root of a clone:   sbatch check_pkg.sh
# Or point it elsewhere:             CROCOTEL_REPO=/path/to/clone sbatch check_pkg.sh
#
# Hoffman2 moved to Rocky 9 + Slurm on 2026-11-05; the Grid Engine version of
# this script (qsub, "#$" directives, -pe shared) no longer runs. -A is
# mandatory because this login belongs to several Slurm accounts. Do NOT add
# --partition: the router picks it from the account.
#SBATCH --job-name=crocotel_check
#SBATCH -A pi_bballiu
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/u/scratch/b/bballiu/crocotel_check_%j.out
#SBATCH --error=/u/scratch/b/bballiu/crocotel_check_%j.err
# NB: variables are NOT expanded in #SBATCH lines, so the log paths above are
# literal. They point at scratch on purpose -- the Grid Engine script used
# "-cwd -j y", which wrote job logs INTO the source tree, and because nothing
# in .Rbuildignore matched them they were built into the v0.3.0 tarball
# (crocotel_check.o14756161 and friends shipped inside the package).

# ~/.profile and the old `mambaload` do not exist on this cluster.
source ~/anaconda3/etc/profile.d/conda.sh
conda activate crocotel
set -eo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# The repo root IS the package (the old $CROCOTEL_REPO/r_package layout was
# archived 2026-10-30). Default to the directory the job was submitted from,
# so any clone checks itself without editing this file.
: "${CROCOTEL_REPO:=${SLURM_SUBMIT_DIR:-$PWD}}"
if [ ! -f "$CROCOTEL_REPO/DESCRIPTION" ]; then
  echo "No DESCRIPTION in CROCOTEL_REPO=$CROCOTEL_REPO" >&2
  echo "Submit from the root of a clone, or set CROCOTEL_REPO to one." >&2
  exit 1
fi
echo "source:  $CROCOTEL_REPO"
echo "version: $(awk '/^Version:/ {print $2}' "$CROCOTEL_REPO/DESCRIPTION")"
R --version | head -1

WORK=/u/scratch/${USER:0:1}/${USER}/pkg_check
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
R CMD build "$CROCOTEL_REPO"
# --no-manual: cluster TeX lacks inconsolata; the manual is checked locally.
R CMD check --as-cran --no-manual crocotel_*.tar.gz
echo "================ 00check.log tail ================"
tail -30 crocotel.Rcheck/00check.log
