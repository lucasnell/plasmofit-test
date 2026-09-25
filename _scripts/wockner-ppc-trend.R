# Option 1: does b_shape ~ 65 actually repair the within-series trend misfit?
#
# The eight-fit panel showed that widening b_shape's prior gains 24.8 elpd and
# shortens cycle_length by 0.47-1.29 h, but sends b_shape from 14.9 to 65.4
# with a posterior sd of 50. Two readings (claude/findings.md, "the eight-fit
# prior panel"):
#
#   (a) the data want a sharp stage window and the default prior suppressed it
#   (b) the model is absorbing structural lack of fit into the one nuisance
#       parameter now free to run
#
# (b) has independent support: the real series rise and fall more steeply than
# the fitted trajectories, median within-series range 2.56 log10 units against
# 1.95-2.21 simulated (claude/findings.md, "Information content"). A sharper
# stage window is exactly what would let the model chase a steeper decline.
#
# So: posterior predictive check on that statistic. If b_shape ~ 65 closes the
# gap to 2.56, the widened fits are describing something the default prior was
# suppressing and reading (a) gains. If the gap survives, the elpd gain is
# being bought somewhere else and reading (b) gains.
#
# WHAT THIS CAN AND CANNOT SHOW. A misfit that persists rules out (a) as a
# complete account. A misfit that closes is consistent with BOTH -- a model
# absorbing lack of fit into a nuisance parameter also fits better. So a
# closed gap does not settle it; that needs the recovery check in
# wockner-schedule-sim.R with SCHEDSIM_TRUTH_FIT=pl_wide_both.
#
# Post-hoc on saved fits, no refitting. Reads four 55 MB fits one at a time.
#   srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#       _scripts/wockner-ppc-trend.R 2>&1 | tee _data/ppc-trend.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

N_DRAW <- 200L
SEED <- 20260925L

## The four fits that span the b_shape contrast, at both log10_total0 settings.
FITS <- c("no_pool", "np_wide_bshape", "np_wide_total0", "np_wide_both")

# ------------------------------------------------------------------------ #
# The real design, built exactly as wockner-fit.R builds it so the group
# codes and series ordering match the saved fits.
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

## Pinned to archer_stan_data()'s current defaults rather than left implicit,
## so a change to the package default cannot silently redefine the baseline.
d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L,
                      mean_log_b_shape = 2, sd_log_b_shape = 0.5,
                      mean_log10_total0 = 1, sd_log10_total0 = 0.25)

ends <- cumsum(d$n_obs)
starts <- ends - d$n_obs + 1L

## The statistic: per-series range of log10(y + 1), then the median over the
## 177 series. Range rather than sd because it is what "rise and fall more
## steeply" means on a series with 4-8 observations, and it is the statistic
## the 2.56 / 1.95-2.21 comparison was made on.
med_range <- function(y) {
    lg <- log10(y + 1)
    median(vapply(seq_len(d$n_ts), \(i) {
        v <- lg[starts[i]:ends[i]]
        max(v) - min(v)
    }, numeric(1)))
}

T_obs <- med_range(d$y)
cat(sprintf("observed median within-series range: %.3f log10 units\n\n", T_obs))

# ------------------------------------------------------------------------ #
# Posterior predictive replicates
#
# Uses the matrix_exp path, as wockner-schedule-sim.R does: it agrees with
# the series solution the fit inverts to 5e-13 (run_check = 1), so this is
# not the code being tested regenerating its own answer.
# ------------------------------------------------------------------------ #

one_draw <- function(p, rng) {
    y_hat <- numeric(d$n_total_obs)
    for (i in seq_len(d$n_ts)) {
        ix <- starts[i]:ends[i]
        y0 <- plasmofit:::generate_starts(p$cycle_length[d$grp_cl[i]], d$n_c,
                                          p$b_shape[d$grp_init[i]],
                                          p$b_offset[d$grp_init[i]],
                                          p$log10_total0[d$grp_init[i]])
        y_hat[ix] <- plasmofit:::mat_exp_series(y0, p$cycle_length[d$grp_cl[i]],
                                                d$n_c, p$R[d$grp_R[i]], d$mu,
                                                d$ts[ix], d$dt_full)
    }
    sd_obs <- p$sd_iRBC[rep(d$grp_sd, d$n_obs)]
    y_rep <- pmax(0, 10^(log10(y_hat + 1) + rnorm(d$n_total_obs, 0, sd_obs)) - 1)
    c(traj = med_range(y_hat), rep = med_range(y_rep))
}

run_fit <- function(cfg) {
    cat("  ", cfg, "... ")
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    pars <- c("b_shape", "b_offset", "log10_total0", "R", "cycle_length",
              "sd_iRBC")
    m <- map(set_names(pars), \(p) as.matrix(f, pars = p))
    n_post <- nrow(m[[1]])
    rm(f); invisible(gc())

    set.seed(SEED)
    ix <- sort(sample.int(n_post, min(N_DRAW, n_post)))
    ## cycle_length under pooled_cl is n_grp_cl identical copies, under
    ## no_pool one per trial; indexing by grp_cl is correct for both.
    out <- map(ix, \(i) {
        p <- map(m, \(x) x[i, ])
        one_draw(p)
    }) |> bind_rows()

    bs <- mean(m$b_shape[ix, ])
    cat(sprintf("done (%d draws, mean b_shape %.1f)\n", length(ix), bs))
    tibble(fit = cfg, b_shape = bs,
           traj_med = median(out$traj), traj_lo = quantile(out$traj, 0.025),
           traj_hi = quantile(out$traj, 0.975),
           rep_med = median(out$rep), rep_lo = quantile(out$rep, 0.025),
           rep_hi = quantile(out$rep, 0.975),
           ppp = mean(out$rep >= T_obs))
}

cat("drawing", N_DRAW, "posterior predictive replicates per fit\n")
res <- map(FITS, run_fit) |> list_rbind()

cat("\nCells: one row per fit. `b_shape` is the posterior mean averaged over\n")
cat("14 grp_init groups. `traj_*` is the median-over-177-series within-series\n")
cat("range of log10(y_hat + 1) for the NOISELESS trajectory, `rep_*` the same\n")
cat("for a full posterior predictive replicate (trajectory + observation\n")
cat("noise), each summarised over 200 draws as median and a 95% interval, in\n")
cat("log10 units. `ppp` is the posterior predictive p-value,\n")
cat("P(T(y_rep) >= T(y_obs)) over draws; ~0.5 is a fit that reproduces the\n")
cat("statistic and near 0 means the model generates series flatter than the\n")
cat(sprintf("real ones. Observed T(y_obs) = %.3f.\n\n", T_obs))

res |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(width = Inf)

write_rds(list(T_obs = T_obs, res = res, n_draw = N_DRAW, seed = SEED),
          "_data/wock-ppc-trend.rds")
cat("\nwrote _data/wock-ppc-trend.rds\n")
