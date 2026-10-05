# Does phase degeneracy at commensurate periods bias an INTEGRAL but not a
# MAXIMUM? That is the signature thread 2 has shown all along.
#
# THE OBSERVATION. 87.8% of within-series intervals are 12 h. At a period
# commensurate with 12 h the observation times land on only a few distinct
# phases: _data/phase-coverage-scan.log puts 7.8 effective phases at the
# simulated truth of 45.012 h and 3.7 at 48 h, with a sharp flat-bottomed dip
# across 47.9-48.1. The estimates drift toward 48.
#
# THE MECHANISM IT WOULD IMPLY. Where phase coverage degenerates, the data
# constrain b_offset weakly, so a wider range of b_offset fits equally well.
# That is parameter-space VOLUME. A posterior mean integrates over volume; a
# maximum does not. findings.md already records exactly that split -- the
# per-trial MLE recovers the truth to -0.195 h while the fitted posterior mean
# sits +1.56 to +1.71 h high -- and it would also explain why widening priors
# never moved the bias (this is not a prior effect) and why the joint surface
# is a comb of modes.
#
# THE TEST, with no MCMC. Fix every parameter at the simulated truth except a
# common phase shift delta applied to all b_offset groups, which is the
# natural "phase of the whole infection" coordinate and the one the sampling
# grid acts on. Then over a grid of cycle_length:
#
#   profile(cl)    = max over delta of the log likelihood
#   integrated(cl) = log of the mean over delta of the likelihood
#
# b_offset's prior in this model is uniform on the circle, so averaging over
# delta with uniform weight IS the correct marginalisation, not a choice.
#
# If profile peaks at the truth and integrated peaks higher -- toward 48 --
# the volume effect is demonstrated directly, and thread 2 has a mechanism.
# If both peak in the same place, it is refuted and the sampling grid is not
# the story.
#
#     srun -N 1 -n 1 -c 4 --mem=32G Rscript --vanilla \
#         _scripts/phase-volume-test.R 2>&1 | tee _data/phase-volume-test.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
Sys.setenv(SLURM_ARRAY_TASK_ID = Sys.getenv("SCHEDSIM_TASK", "1"))
src <- readLines("_scripts/wockner-schedule-sim.R")
cut <- grep("^# Fit$", src) - 2L
eval(parse(text = paste(src[seq_len(cut)], collapse = "\n")))   # gives d, truth, y_sim

suppressPackageStartupMessages({ library(tidyverse) })

## Grid size matters and was sized AFTER launching the first time, which is
## the wrong order: each evaluation builds a 2*n_c = 192 square matrix and
## exponentiates it once per series, so the cost is periods x phases x 177
## matrix exponentials. The first grid, 0.05 h by 1 degree, is ~72,000 points
## and runs for most of a day. 0.25 h by 5 degrees is ~2,400 and resolves the
## only question asked here: whether two peaks sit near 45 or near 48.
CL_GRID <- seq(42, 50, by = 0.25)
N_DELTA <- 72                       # phase resolution, 5 degrees

ends <- cumsum(d$n_obs); starts <- ends - d$n_obs + 1L
y_obs <- log10(y_sim + 1)
sd_obs <- truth$sd_iRBC[rep(d$grp_sd, d$n_obs)]

## log likelihood at a given cycle_length and a common phase shift
ll <- function(cl, delta) {
    bo <- (truth$b_offset + delta) %% 1
    yh <- numeric(d$n_total_obs)
    for (i in seq_len(d$n_ts)) {
        ix <- starts[i]:ends[i]
        y0 <- plasmofit:::generate_starts(cl, d$n_c, truth$b_shape[d$grp_init[i]],
                                          bo[d$grp_init[i]],
                                          truth$log10_total0[d$grp_init[i]])
        yh[ix] <- plasmofit:::mat_exp_series(y0, cl, d$n_c, truth$R[d$grp_R[i]],
                                             d$mu, d$ts[ix], d$dt_full)
    }
    sum(dnorm(y_obs, log10(yh + 1), sd_obs, log = TRUE))
}

deltas <- seq(0, 1, length.out = N_DELTA + 1)[-(N_DELTA + 1)]
cat("=== scanning", length(CL_GRID), "periods x", N_DELTA, "phases ===\n")
res <- map(CL_GRID, \(cl) {
    v <- vapply(deltas, \(dd) ll(cl, dd), numeric(1))
    m <- max(v)
    tibble(cl = cl, profile = m,
           integrated = m + log(mean(exp(v - m))),   # log mean exp, stable
           delta_hat = deltas[which.max(v)])
}) |> list_rbind()

cat("\nCells: `profile` is the maximum of the log likelihood over the common
phase shift at that cycle_length; `integrated` is the log of its mean over
phase with uniform weight, which is the correct marginalisation because
b_offset's prior here is uniform on the circle. Both in log units, relative
to their own maximum so only shape matters. `volume` is integrated minus
profile -- how much phase-space is compatible with the data at that period,
in log units, higher meaning more. Every parameter but the phase shift is
held at the simulated truth.\n\n")

out <- res |> mutate(profile_rel = profile - max(profile),
                     integrated_rel = integrated - max(integrated),
                     volume = integrated - profile)
out |> filter(cl %in% c(42, 43, 44, 45, 46, 47, 47.5, 48, 48.5, 49, 50)) |>
    transmute(cl, profile_rel = round(profile_rel, 2),
              integrated_rel = round(integrated_rel, 2),
              volume = round(volume, 2)) |> print(n = Inf)

cat(sprintf("\n  true cycle_length          %.3f h\n", truth$cycle_length))
cat(sprintf("  PROFILE likelihood peaks at %.2f h  (bias %+.2f h)\n",
            out$cl[which.max(out$profile)],
            out$cl[which.max(out$profile)] - truth$cycle_length))
cat(sprintf("  INTEGRATED likelihood peaks at %.2f h  (bias %+.2f h)\n",
            out$cl[which.max(out$integrated)],
            out$cl[which.max(out$integrated)] - truth$cycle_length))
cat(sprintf("  phase volume at 48 h minus at the truth: %+.2f log units\n",
            out$volume[which.min(abs(out$cl - 48))] -
            out$volume[which.min(abs(out$cl - truth$cycle_length))]))
write_rds(out, sprintf("_data/wock-phase-volume-%s.rds", cfg_name))
