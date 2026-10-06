#!/bin/bash -l

# Thread 5 on the PRODUCTION config: does the n_c effect survive a pinned
# b_shape?
#
# 29633 (entries 37-38) ran with b_shape ESTIMATED, which is the right setting
# for asking whether b_shape absorbs the desynchronisation rate but is NOT the
# configuration this project has settled on. It found:
#   - b_shape falls monotonically, 64.5 -> 61.6 -> 25.9 across n_c 96/192/384
#   - cycle_length moves -3.20 h, 44.25 -> 42.07 -> 41.04
#   - n_c = 96 is 71.8 elpd worse than 384 (se 13.2); 192 vs 384 is -3.2 (se
#     5.0), so the ladder plateaus at 192
#   - max_rel_diff at n_c = 96 is ~1e-12, so the arithmetic is faithful and the
#     effect is STRUCTURAL, not numerical (_scripts/nc-numerical-check.R)
#
# Entries 39-40 pin b_shape at 400 -- the setting that gains +28 to +37 elpd on
# real data -- against the n_c = 96 baseline np_bs400_data (entry 22). This
# decides whether the headline cycle-length number is affected.
#
# PREDICTION, recorded before the fits run and repeated in wockner-fit.R:
#   - effect LARGELY GONE    -> the sensitivity was b_shape absorbing a wrong
#                               decay rate; the production estimate stands
#   - effect PERSISTS        -> n_c is misspecified independently of b_shape,
#                               the production estimate is off by ~3 h, and
#                               every cycle-length number needs requalifying
#   - effect LARGER          -> b_shape was partly compensating, and pinning it
#                               exposes more of the error
#
# Cost, measured on 29633 rather than extrapolated: 5h47 at n_c = 192 and
# 15h59 at 384. The nc-sizing.R probe underestimated because it measured cost
# per leapfrog (1.81x) and the leapfrog count ALSO rose, 223 -> 390.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-fit-nc-bs400.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=39-40
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=wock-ncbs
#SBATCH --chdir=/home2/lan68/plasmofit/plasmofit-test/_data
#SBATCH --output=wock-ncbs-%a.out
#SBATCH --error=wock-ncbs-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla ../_scripts/wockner-fit.R
