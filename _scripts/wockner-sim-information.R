# Why do the simulated data identify the nuisance parameters less well than
# the real data, at the same observation times and the same sample size?
#
# Open thread 3. The simulation returns b_shape ~25-40% low and log10_total0
# ~100% high whatever the truth is built from, while the real fit moves both
# far from their priors (b_shape 7.39 -> 18.9, log10_total0 1 -> 0.312). The
# lead is that simulated log10(y + 1) has sd 0.85-0.89 against a real 1.01,
# so the simulated observations vary less than the real ones.
#
# This script is post-hoc: it reads saved fits and regenerates the simulated
# datasets from their seeds. No new sampling.
#
# The quantity that matters is not the spread of the data but the spread of
# the SIGNAL relative to the noise. Information about a periodic parameter
# comes from the oscillation, not from the overall trend, so each series is
# detrended on time and what is left is the oscillatory component:
#
#   wiggle_i  = sd over series i of (log10(y_hat + 1) detrended on t)
#   noise_i   = sd_iRBC for that series' grp_sd group
#   snr_i     = wiggle_i / noise_i
#
# snr^2 summed over observations is the usual Fisher-information scaling for
# a signal in additive noise, so it is reported that way.
#
# Runs where the fits are; loads two full fits.
#   srun -N 1 -n 1 -c 4 --mem=32G Rscript --vanilla \
#       _scripts/wockner-sim-information.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

pm <- function(f, par) unname(colMeans(as.matrix(f, pars = par)))


# ------------------------------------------------------------------------ #
# The design, identical for the real data and every simulated replicate
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L)

ends <- cumsum(d$n_obs)
starts <- ends - d$n_obs + 1L

## Noiseless trajectories at a given parameter set. cycle_length may be a
## scalar (the simulation gives every trial the same one) or one per grp_cl
## (the real fit); recycled either way.
traj <- function(p) {
    cl <- if (length(p$cycle_length) == 1L) rep(p$cycle_length, d$n_grp_cl)
          else p$cycle_length
    out <- numeric(d$n_total_obs)
    for (i in seq_len(d$n_ts)) {
        ix <- starts[i]:ends[i]
        y0 <- plasmofit:::generate_starts(cl[d$grp_cl[i]], d$n_c,
                                          p$b_shape[d$grp_init[i]],
                                          p$b_offset[d$grp_init[i]],
                                          p$log10_total0[d$grp_init[i]])
        out[ix] <- plasmofit:::mat_exp_series(y0, cl[d$grp_cl[i]], d$n_c,
                                              p$R[d$grp_R[i]], d$mu,
                                              d$ts[ix], d$dt_full)
    }
    out
}

## Residual sd after removing a within-series linear trend in time. Series
## with fewer than 3 observations carry no residual after a 2-parameter fit
## and are dropped rather than scored 0.
detrended_sd <- function(v) {
    map_dbl(seq_len(d$n_ts), \(i) {
        ix <- starts[i]:ends[i]
        if (length(ix) < 3L) return(NA_real_)
        stats::sd(stats::residuals(stats::lm(v[ix] ~ d$ts[ix])))
    })
}

series_sd_iRBC <- function(sdv) sdv[d$grp_sd]

## Everything reported for one dataset: the observed values, the noiseless
## trajectory they came from, and the per-group observation-error sds.
describe <- function(label, y_obs, y_hat, sd_iRBC) {
    lg  <- log10(y_obs + 1)
    lgh <- log10(y_hat + 1)
    series_of <- rep(seq_len(d$n_ts), d$n_obs)
    series_mean <- tapply(lg, series_of, mean)
    within_sd <- tapply(lg, series_of, stats::sd)

    wig <- detrended_sd(lgh)          # signal: oscillation after detrending
    obs_wig <- detrended_sd(lg)       # signal + noise, for reference
    nz <- series_sd_iRBC(sd_iRBC)
    snr <- wig / nz
    n_per <- d$n_obs

    keep <- is.finite(snr)
    tibble(
        dataset = label,
        sd_total = stats::sd(lg),
        sd_between_series = stats::sd(series_mean),
        sd_within_series = mean(within_sd, na.rm = TRUE),
        range_med = median(tapply(lg, series_of, \(x) diff(range(x)))),
        wiggle_med = median(wig, na.rm = TRUE),
        obs_wiggle_med = median(obs_wig, na.rm = TRUE),
        noise_med = median(nz),
        snr_med = median(snr, na.rm = TRUE),
        info = sum((snr[keep]^2) * n_per[keep])
    )
}


# ------------------------------------------------------------------------ #
# Real data, at the real no_pool fit's posterior means
# ------------------------------------------------------------------------ #

f_np <- read_rds("_data/wock-fit-no_pool.rds")
p_real <- list(cycle_length = pm(f_np, "cycle_length"),
               b_shape = pm(f_np, "b_shape"),
               b_offset = pm(f_np, "b_offset"),
               log10_total0 = pm(f_np, "log10_total0"),
               R = pm(f_np, "R"),
               sd_iRBC = pm(f_np, "sd_iRBC"))
rm(f_np); invisible(gc())

rows <- list(describe("real", d$y, traj(p_real), p_real$sd_iRBC))


# ------------------------------------------------------------------------ #
# Simulated replicates, regenerated from their seeds
#
# The truth is rebuilt the way wockner-schedule-sim.R builds it: the
# pooled_cl fit's posterior mean vector, or a single posterior draw. The
# noise draw is the replicate's own seed, so these are the same datasets the
# fits saw, not new ones.
# ------------------------------------------------------------------------ #

f_pl <- read_rds("_data/wock-fit-pooled_cl.rds")
pars <- c("b_shape", "b_offset", "log10_total0", "R", "sd_iRBC", "cycle_length")
mean_truth <- set_names(map(pars, \(p) pm(f_pl, p)), pars)
mean_truth$cycle_length <- mean_truth$cycle_length[1]
draw_of <- function(i) {
    o <- set_names(map(pars, \(p) unname(as.matrix(f_pl, pars = p)[i, ])), pars)
    o$cycle_length <- o$cycle_length[1]
    o
}

REP_SEEDS <- c(742160183L, 1908445027L, 355721694L)
sims <- tribble(
    ~label,     ~rep, ~draw,
    "sim-rep1",    1L,   NA_integer_,
    "sim-rep2",    2L,   NA_integer_,
    "sim-rep3",    3L,   NA_integer_,
    "sim-draw500", 2L,  500L,
    "sim-draw1500",2L, 1500L,
    "sim-draw2500",2L, 2500L
)

for (k in seq_len(nrow(sims))) {
    tr <- if (is.na(sims$draw[k])) mean_truth else draw_of(sims$draw[k])
    y_hat <- traj(tr)
    set.seed(REP_SEEDS[sims$rep[k]])
    sd_obs <- tr$sd_iRBC[rep(d$grp_sd, d$n_obs)]
    y_sim <- pmax(0, 10^(log10(y_hat + 1) + rnorm(d$n_total_obs, 0, sd_obs)) - 1)
    rows <- c(rows, list(describe(sims$label[k], y_sim, y_hat, tr$sd_iRBC)))
}
rm(f_pl); invisible(gc())

out <- list_rbind(rows)


# ------------------------------------------------------------------------ #
# Report
# ------------------------------------------------------------------------ #

cat("=== spread of log10(y + 1) ===\n")
cat("  sd_total decomposes into variation BETWEEN series and WITHIN them.\n")
out |> select(dataset, sd_total, sd_between_series, sd_within_series,
              range_med) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat("\n=== signal against noise ===\n")
cat("  wiggle_med: median over series of the detrended sd of the NOISELESS\n")
cat("    trajectory -- the oscillation a periodic parameter is read from.\n")
cat("  noise_med:  median observation-error sd on the same scale.\n")
cat("  info:       sum over observations of (wiggle / noise)^2.\n")
out |> select(dataset, wiggle_med, obs_wiggle_med, noise_med, snr_med, info) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

real <- out |> filter(dataset == "real")
cat("\n=== simulated relative to real ===\n")
out |> filter(dataset != "real") |>
    transmute(dataset,
              sd_total = sd_total / real$sd_total,
              between = sd_between_series / real$sd_between_series,
              within = sd_within_series / real$sd_within_series,
              wiggle = wiggle_med / real$wiggle_med,
              noise = noise_med / real$noise_med,
              info = info / real$info) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat("\n  A ratio near 1 on `info` means the simulated design carries as much\n",
    "  information as the real one and the mis-recovery needs another\n",
    "  explanation. A ratio well under 1 locates it: read `wiggle` and\n",
    "  `noise` to see whether the signal shrank or the noise grew.\n", sep = "")

write_rds(out, "_data/wock-sim-information.rds")
cat("\nwrote _data/wock-sim-information.rds\n")
