#!/bin/bash -l

# Thread 15, the decay-law test: does synchrony decay as sqrt(cycles) or
# linearly in cycles?
#
# The Erlang chain gives sd = cycle_length * sqrt(k / n_c), a DIFFUSIVE sqrt
# law, and a 1-D IPM with a Gaussian development kernel gives the same law.
# So switching to an IPM would decouple the dispersion rate from the
# numerical mesh -- the defect diagnosed in findings.md -- WITHOUT testing
# whether the functional form is right. The competing form is
# between-parasite heterogeneity in cycle duration, which spreads LINEARLY in
# cycles. This test decides which, and therefore decides what an IPM's kernel
# should be, so it is ordered strictly before any rewrite.
#
# Maximum likelihood, not Bayes: a screen whose only job is to say whether
# the Stan work is worth doing. It exploits the fact that b_shape, b_offset
# and log10_total0 are per grp_init while R and cycle_length are per trial,
# so within a grp_init unit every series shares every parameter and the unit
# is one trajectory with its observations as replicates. 14 units, one task
# each, independent.
#
# Cost is roughly CUBIC in n_c (0.0148 s per trajectory at 96, 0.153 at 192,
# 1.18 at 384, 9.21 at 768, measured), and model B multiplies it by the 7
# quadrature nodes, which is why this is an array rather than one job.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/decay-law-test.sh
#   Rscript --vanilla _scripts/decay-law-read.R

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-14
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --job-name=decay-law
#SBATCH --output=_data/decay-law-%a.out
#SBATCH --error=_data/decay-law-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

Rscript --vanilla _scripts/decay-law-test.R
