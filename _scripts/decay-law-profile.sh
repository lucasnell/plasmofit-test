#!/bin/bash -l

# Thread 15, the prerequisite named in claude/ipm-decision.md: is the
# dispersion MAGNITUDE identified at all, with the initial-spread nuisance
# free?
#
# WHY THIS IS THE RIGHT INSTRUMENT. In the chain, transit coefficient of
# variation is 1/sqrt(n_c), so n_c IS the dispersion magnitude, exactly and
# monotonically. b_shape is free at every rung, so profiling the maximised
# log-likelihood over n_c is a PROFILE LIKELIHOOD in the dispersion with the
# initial-spread nuisance concentrated out -- which is the identifiability
# question, asked with code that is already built and validated. Staying
# inside the sqrt family costs nothing: 29684 showed the form is not learnable
# from these data, only the magnitude.
#
# WHAT IT ADDS. n_c = 128, 256, 512, filling the gaps between the 96/192/384
# already run (29684) and the 768 running now (30524). Seven rungs spanning a
# transit CV of 0.102 down to 0.036 is enough to show CURVATURE rather than
# just an ordering, and curvature is what an interval needs.
#
# WHY IT MATTERS FOR THE IPM. If the pooled profile has an interior maximum
# with a finite 2-log-likelihood interval, the dispersion rate is a MEASURED
# QUANTITY and an IPM would deliver a parameter with an interval -- the
# strongest argument for building one. If the profile runs monotonically to
# the fine end, "free dispersion" means "dispersion goes to zero", and an IPM
# buys a boundary estimate rather than a rate. That is worth knowing BEFORE
# discarding the Erlang-window series.
#
# ALSO RECORDED: the fitted b_shape at every rung, so the b_shape/dispersion
# trade-off -- the actual confound -- can be read directly.
#
# COST. Cubic in n_c: (128^3 + 256^3 + 512^3) / 384^3 = 2.70 times model A at
# 384, which took 1777 s per unit in 29684, so ~1.3 h per unit and 1-3 h
# across units, which differ about twofold.
#
# RESOURCES. 14 x 1 CPU x 8G = 14 CPUs and 114,688 MB. Running alongside
# 30524 (14 CPUs, 229,376 MB) and 30525 (4 CPUs, 24,576 MB) that totals 32
# CPUs and 368,640 MB, which is 25% of the 128-CPU budget and 71% of the
# 515,670 MB budget -- inside the half-node rule in CLAUDE.md, so no
# --array=1-N%M throttle is needed. Check squeue before adding a fourth job.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/decay-law-profile.sh
#   Rscript --vanilla _scripts/nc-dispersion-profile-read.R | tee _data/nc-dispersion-profile-$(date +%F).txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=decayprof
#SBATCH --output=_data/decay-law-prof-%a.out
#SBATCH --error=_data/decay-law-prof-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

export DECAY_MODE=profile
Rscript --vanilla _scripts/decay-law-test.R
