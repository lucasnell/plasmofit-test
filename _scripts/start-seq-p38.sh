#!/bin/bash -l

# Rerun of _scripts/start-seq.R task 38 (entry 59, trial SJ733IBSMCS), which
# hit the 24 h wall at 14 of 16 p3 points having saved nothing. One array task
# per p3 grid point (SSQ_P3ONLY), so each is ~2 h instead of ~27 h; merge mode
# stitches the per-point files automatically. Everything else is identical to
# the main run, except that a single-point task has no warm start from the
# previous grid point.
#
# FOOTPRINT. 2 CPUs and 2G per task, 16 tasks at once: 32 CPUs and 32,768 MB,
# well inside the 128 CPU / 515,670 MB half-node budget.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/start-seq-p38.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-16
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --job-name=start-seq-p38
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/start-seq/slurm-p38-%a.out
#SBATCH --error=_data/start-seq/slurm-p38-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SSQ_MODE=run SSQ_P3ONLY=${SLURM_ARRAY_TASK_ID} SLURM_ARRAY_TASK_ID=38 \
    Rscript --vanilla _scripts/start-seq.R
