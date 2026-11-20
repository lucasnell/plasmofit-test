#!/bin/bash -l

# Profile likelihood of the sequestration midpoint p3 under the starting
# distributions of entries 54, 57 and 59; design and reading rule in
# _scripts/start-seq.R. 39 tasks = 3 models x 13 trials, each one
# single-threaded R process, no Stan sampling. 2 CPUs per task so each has a
# whole physical core (hyperthread sibling idle).
#
# FOOTPRINT. 2 CPUs and 2G per task, at most 32 at once: 64 CPUs and 65,536 MB,
# inside the 128 CPU / 515,670 MB half-node budget
# (M = floor(min(128/2, 515670/2048)) = 64; %32 leaves half of it).
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/start-seq.sh
#   SSQ_MODE=merge Rscript --vanilla _scripts/start-seq.R > _data/start-seq-<date>.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-39%32
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=start-seq
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test
#SBATCH --output=_data/start-seq/slurm-%a.out
#SBATCH --error=_data/start-seq/slurm-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
SSQ_MODE=run Rscript --vanilla _scripts/start-seq.R
