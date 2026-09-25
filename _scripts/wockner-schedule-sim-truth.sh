#!/bin/bash -l

# Option 2: does this design recover b_shape = 65, or send it back to ~15?
#
# The eight-fit panel put b_shape at 65.4 (posterior sd 50) once its prior was
# widened, and gained 24.8 elpd doing it. Two readings (claude/findings.md,
# "the eight-fit prior panel"): the data want a sharp stage window, or the
# model is absorbing lack of fit into the one nuisance parameter free to run.
#
# This simulates from `pl_wide_both` -- the pooled fit WITH both priors
# widened, so its truth has b_shape ~ 65 -- and fits under those same widened
# priors. If the fit returns ~15, the real-data 65 is not a measurement and
# the absorption reading gains. If it returns ~65, the design can see a sharp
# window when one is there, and the real-data 65 has to be taken seriously.
#
# Arm `wide_nuis` matches pl_wide_both's priors exactly (sd_log_b_shape 1.5,
# sd_log10_total0 1), so the fit is not being handed a prior the truth was
# not generated under. Tasks 31-33 are wide_nuis replicates 1-3; the truth
# override renames the outputs to `wide_nuis-rep<n>-truthpl_wide_both`, so
# nothing collides with the wide_nuis fits already on disk.
#
# NOTE the truth's cycle_length is 43.7 h here, not the 45.012 h every other
# simulation in this project uses. Read bias against what the job prints.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/wockner-schedule-sim-truth.sh

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=31-33
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=wock-schedsim-truth
#SBATCH --output=_data/wock-schedsim-truth-%a.out
#SBATCH --error=_data/wock-schedsim-truth-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export SCHEDSIM_TRUTH_FIT=pl_wide_both

Rscript --vanilla _scripts/wockner-schedule-sim.R
