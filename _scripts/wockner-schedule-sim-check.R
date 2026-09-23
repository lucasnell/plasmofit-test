# Is the schedule-bias simulation's upward cycle_length bias real, or does the
# simulator simply disagree with the likelihood?
#
# wockner-schedule-sim.R generates trajectories with mat_exp_series and then
# fits a model whose likelihood evaluates them with ew_poly_series, a different
# solver. Every replicate came back biased +1.2 to +2.2 h above the known
# truth. Before reading that as a property of the observation schedules, rule
# out the mundane explanation: that the two solvers do not agree, so the fit is
# inverting data no model generated.
#
# archer_fit.stan already carries the cross-check. With run_check = 1 its
# generated quantities computes max_rel_diff, the largest relative gap between
# ew_poly_series and full_mat_exp_series over every series and observation
# time, evaluated at each draw's parameters. mat_exp_series and
# full_mat_exp_series are the same solver read at different indices (verified
# at the ts/dt_full + 1 offset), so a max_rel_diff near machine precision means
# the simulator and the likelihood agree and the bias is not an artifact.
#
# Sampling quality is irrelevant here -- the check is a property of the solvers
# at whatever parameters the sampler visits -- so this runs short and cheap.
#
#   srun -N 1 -n 1 -c 2 --mem=16G --time=01:00:00 \
#       Rscript --vanilla _scripts/wockner-schedule-sim-check.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

## Replicate 1 of the default arm, rebuilt exactly: same seed, same truth.
SEED_NOISE <- 742160183L
SEED_FIT <- 538065874L

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L,
                      run_check = 1L)

post_mean <- function(f, par) unname(colMeans(as.matrix(f, pars = par)))

f_pl <- read_rds("_data/wock-fit-pooled_cl.rds")
truth <- list(
    b_shape      = post_mean(f_pl, "b_shape"),
    b_offset     = post_mean(f_pl, "b_offset"),
    log10_total0 = post_mean(f_pl, "log10_total0"),
    R            = post_mean(f_pl, "R"),
    sd_iRBC      = post_mean(f_pl, "sd_iRBC")
)
truth$cycle_length <- post_mean(f_pl, "cycle_length")[1]
rm(f_pl); invisible(gc())

ends <- cumsum(d$n_obs)
starts <- ends - d$n_obs + 1L

y_hat <- numeric(d$n_total_obs)
for (i in seq_len(d$n_ts)) {
    ix <- starts[i]:ends[i]
    y0 <- plasmofit:::generate_starts(truth$cycle_length, d$n_c,
                                      truth$b_shape[d$grp_init[i]],
                                      truth$b_offset[d$grp_init[i]],
                                      truth$log10_total0[d$grp_init[i]])
    y_hat[ix] <- plasmofit:::mat_exp_series(y0, truth$cycle_length, d$n_c,
                                            truth$R[d$grp_R[i]], d$mu,
                                            d$ts[ix], d$dt_full)
}

set.seed(SEED_NOISE)
sd_obs <- truth$sd_iRBC[rep(d$grp_sd, d$n_obs)]
d$y <- pmax(0, 10^(log10(y_hat + 1) + rnorm(d$n_total_obs, 0, sd_obs)) - 1)

f <- archer_fit(d, model = "no_pool", chains = 2L, iter = 150L, warmup = 75L,
                seed = SEED_FIT, threads_per_chain = 1L)

mrd <- as.numeric(rstan::extract(f, "max_rel_diff")[[1]])

cat("\n=== ew_poly_series vs full_mat_exp_series on the simulated data ===\n")
cat(sprintf("  max_rel_diff over %d draws: max %.3e, median %.3e\n",
            length(mrd), max(mrd), median(mrd)))
cat(sprintf("  draws exceeding 1e-6: %d of %d\n", sum(mrd > 1e-6), length(mrd)))

if (max(mrd) < 1e-6) {
    cat("\n  The two solvers agree to better than 1e-6 everywhere. The\n",
        "  simulator and the likelihood share a forward model, so the\n",
        "  cycle_length bias is a property of the design, not a solver\n",
        "  mismatch.\n", sep = "")
} else {
    cat("\n  WARNING: the solvers disagree. The bias reported by\n",
        "  wockner-schedule-sim-analyze.R may be an artifact of generating\n",
        "  with one solver and fitting with another. Do not interpret it\n",
        "  until this is resolved.\n", sep = "")
}
