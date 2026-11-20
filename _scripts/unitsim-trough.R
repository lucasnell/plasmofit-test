## Trough-depth posterior predictive check for the single-unit fits
## (2026-11-02, Lucas). elpd_loo sums over every observation and is dominated
## by the overall level, which all candidates get about right; what the fixed
## (b_shape, n_c) candidates actually differ in is the DEPTH OF THE TROUGHS in
## circulating parasites, since sequestration hides mature stages and a
## narrower stage spread cuts a deeper trough. This checks that feature alone.
##
## STATISTIC, fixed before any fit was scored. On y_trans = log10(y + 1) over
##   the unit's observations, remove an OLS straight line in time (pooling all
##   series of the unit, as they share one trajectory), then
##     T_range = max(resid) - min(resid)   the peak-to-trough range
##     T_sd    = sd(resid)
##   computed identically on the observed data and on each posterior
##   predictive replicate. p = P(T_rep >= T_obs) over draws. A candidate whose
##   stage spread is too NARROW oscillates too deeply, so T_rep > T_obs and p
##   is near 1; too wide gives p near 0; p near 0.5 is calibrated. The
##   selection statistic is therefore |p - 0.5|, smaller being better.
##
## HOW THE PREDICTIVE MEAN IS RECOVERED, without re-running any trajectory.
##   The model's log_lik is normal_lpdf(y_trans[k] | mu_k, sd_iRBC), so
##     |y_trans[k] - mu_k| = sd * sqrt(-2 (log_lik + log(sd) + log(2 pi)/2)).
##   That leaves a sign, but within one grp_init unit every series shares one
##   trajectory, so all observations at the same time share one mu; with 4-6
##   observations per time the consistent root is unique. Checked against a
##   direct generate_starts + mat_exp_series trajectory in check mode.
##   Replicates then use the same draw of sd_iRBC and the same truncation the
##   data were simulated with: y_rep_trans = max(0, mu + eps).
##
## MODES (UNITSIM_TROUGH):
##   check    recovered mu against the direct trajectory, 3 draws of one fit.
##   compute  one trough-<fit>.rds per matching fit (skips existing).
##            UNITSIM_TROUGH_MATCH, a regex, restricts which fits (default:
##            the fixed-pin grid and its reference, i.e. gen 384 fix fits
##            under the truth-centred prior).
##   read     selection by |p - 0.5| against elpd, the least-biased candidate,
##            and a random pick; keyed as _scripts/unitsim-pins-read.R.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   UNITSIM_TROUGH=check   Rscript --vanilla _scripts/unitsim-trough.R
##   UNITSIM_TROUGH=compute Rscript --vanilla _scripts/unitsim-trough.R
##   UNITSIM_TROUGH=read    Rscript --vanilla _scripts/unitsim-trough.R > _data/unitsim-read-trough-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({
    library(plasmofit); library(rstan); library(dplyr); library(tidyr)
})
options(width = 165)
D <- "_data/unitsim/study/"
MODE <- Sys.getenv("UNITSIM_TROUGH", "read")
MATCH <- Sys.getenv("UNITSIM_TROUGH_MATCH", "-g384-fix-rep[0-9]+(-pin[a-z0-9]+)?(-f(288|512))?-ptruth$")
SEED <- 20261102L
mod <- plasmofit:::stanmodels[["archer_fit_pooled_both"]]
## Data lists saved before plasmofit's start_dist build (live 2026-11-18)
## lack its three fields, which the Stan programs now require; give them
## the old behaviour (claude/gotchas.md).
with_start_fields <- function(d) {
    if (is.null(d$start_dist)) {
        d$start_dist <- 0L; d$start_window <- 1; d$start_shape_cap <- 3
    }
    ## start_below/start_delay (plasmofit build of 2026-11-20); off values,
    ## ignored by builds that do not declare them
    if (is.null(d$start_below)) {
        d$start_below <- 0; d$start_delay <- 0L; d$delay_max <- 0
    }
    d
}
trough_path <- function(nm) sprintf("%strough-%s.rds", D, nm)

## resid of y_trans after removing an OLS line in time; the design matrix is
## the same for observed and replicated data, so build it once per fit.
detrend_mat <- function(ts) { X <- cbind(1, ts); diag(length(ts)) - X %*% solve(crossprod(X), t(X)) }

## mu per draw and distinct time, from log_lik (see header).
recover_mu <- function(ll, sd, ytr, ts) {
    u <- sort(unique(ts))
    q <- -2 * (ll + log(sd) + 0.5 * log(2 * pi)); q[q < 0] <- 0
    del <- sd * sqrt(q)                                   # draws x obs
    mu <- matrix(NA_real_, nrow(ll), length(u))
    for (j in seq_along(u)) {
        k <- which(ts == u[j]); stopifnot(length(k) >= 2L)
        lo <- ytr[k[1]] - del[, k[1]]; hi <- ytr[k[1]] + del[, k[1]]
        err <- function(cand) {                           # discrepancy of a root
            e <- 0
            for (i in k[-1]) e <- e + pmin(abs(cand - (ytr[i] - del[, i])),
                                           abs(cand - (ytr[i] + del[, i])))
            e
        }
        e_lo <- err(lo); e_hi <- err(hi)
        mu[, j] <- ifelse(e_lo <= e_hi, lo, hi)
        bad <- pmin(e_lo, e_hi) > 1e-4
        if (any(bad)) warning(sprintf("time %g: %d draws with ambiguous root", u[j], sum(bad)))
    }
    list(mu = mu, u = u)
}

stat_of <- function(y, Q) { r <- as.vector(Q %*% y); c(range = max(r) - min(r), sd = sd(r)) }

trough_of <- function(nm) {
    f <- readRDS(sprintf("%sfit-%s.rds", D, nm))
    d <- readRDS(sprintf("%schain-%s-c1.rds", D, nm))$d
    d$calc_log_lik <- 1L
    d <- with_start_fields(d)
    dr <- as.matrix(f)
    g <- rstan::gqs(mod, data = d, draws = dr)
    ll <- rstan::extract(g, "log_lik")$log_lik
    sd <- dr[, "sd_iRBC[1]"]
    ytr <- log10(d$y + 1)
    rm <- recover_mu(ll, sd, ytr, d$ts)
    Q <- detrend_mat(d$ts)
    obs <- stat_of(ytr, Q)
    idx <- match(d$ts, rm$u)
    set.seed(SEED)
    n <- length(ytr); nd <- nrow(ll)
    rep_stat <- matrix(NA_real_, nd, 2L)
    for (i in seq_len(nd)) {
        yrep <- pmax(0, rm$mu[i, idx] + rnorm(n, 0, sd[i]))
        rep_stat[i, ] <- stat_of(yrep, Q)
    }
    list(name = nm, n_obs = n, draws = nd,
         T_range_obs = obs[["range"]], T_sd_obs = obs[["sd"]],
         T_range_rep_mean = mean(rep_stat[, 1]), T_sd_rep_mean = mean(rep_stat[, 2]),
         p_range = mean(rep_stat[, 1] >= obs[["range"]]),
         p_sd = mean(rep_stat[, 2] >= obs[["sd"]]))
}

if (MODE == "check") {
    nm <- "w192-g384-fix-rep1-ptruth"
    f <- readRDS(sprintf("%sfit-%s.rds", D, nm))
    d <- readRDS(sprintf("%schain-%s-c1.rds", D, nm))$d; d$calc_log_lik <- 1L
    d <- with_start_fields(d)
    dr <- as.matrix(f)[1:3, , drop = FALSE]
    g <- rstan::gqs(mod, data = d, draws = dr)
    ll <- rstan::extract(g, "log_lik")$log_lik
    rm <- recover_mu(ll, dr[, "sd_iRBC[1]"], log10(d$y + 1), d$ts)
    cat("Recovered mu (from log_lik) against a direct generate_starts + mat_exp_series\n",
        "trajectory, log10(y + 1) at the unit's distinct observation times.\n\n", sep = "")
    for (i in 1:3) {
        y0 <- plasmofit:::generate_starts(dr[i, "cycle_length[1]"], d$n_c, d$b_shape_data[1],
                                          dr[i, "b_offset[1]"], dr[i, "log10_total0[1]"])
        yh <- log10(plasmofit:::mat_exp_series(y0, dr[i, "cycle_length[1]"], d$n_c,
                                               dr[i, "R[1]"], d$mu, rm$u, d$dt_full) + 1)
        cat(sprintf(" draw %d  max |recovered - direct| %.3e\n", i, max(abs(rm$mu[i, ] - yh))))
        stopifnot(max(abs(rm$mu[i, ] - yh)) < 1e-6)
    }
    x <- trough_of(nm)
    cat(sprintf("\n%s: T_range observed %.3f, replicate mean %.3f, p %.3f; T_sd observed %.3f, replicate mean %.3f, p %.3f\n",
                nm, x$T_range_obs, x$T_range_rep_mean, x$p_range, x$T_sd_obs, x$T_sd_rep_mean, x$p_sd))
    quit(save = "no")
}

if (MODE == "compute") {
    nms <- sub("^.*sum-(.*)\\.rds$", "\\1", Sys.glob(paste0(D, "sum-*.rds")))
    nms <- nms[grepl(MATCH, nms)]
    cat(length(nms), "fits match\n")
    n <- 0L
    for (nm in nms) {
        if (file.exists(trough_path(nm))) next
        saveRDS(trough_of(nm), trough_path(nm))
        n <- n + 1L
        if (n %% 50L == 0L) cat(format(Sys.time()), n, "done\n")
    }
    cat("computed", n, "this run;", length(Sys.glob(paste0(D, "trough-*.rds"))), "in all\n")
    quit(save = "no")
}

if (MODE != "read") stop("unknown UNITSIM_TROUGH: ", MODE)

## ---- read ------------------------------------------------------------------
s <- bind_rows(lapply(Sys.glob(paste0(D, "sum-*.rds")), readRDS))
for (v in c("prior", "truth_set", "pin")) if (!v %in% names(s)) s[[v]] <- NA_character_
if (!"fit_nc" %in% names(s)) s$fit_nc <- NA_integer_
s <- s |> mutate(prior = coalesce(prior, "prod"), truth_set = coalesce(truth_set, "main"),
                 pin = coalesce(pin, "truth"), fit_nc = coalesce(fit_nc, 384L),
                 cell = paste0(pin, "-f", fit_nc))
tf <- Sys.glob(paste0(D, "trough-*.rds"))
tg <- lapply(tf, readRDS); names(tg) <- vapply(tg, `[[`, character(1), "name")
s <- filter(s, name %in% names(tg)) |>
    mutate(p_range = vapply(tg[name], `[[`, numeric(1), "p_range"),
           p_sd = vapply(tg[name], `[[`, numeric(1), "p_sd"),
           T_obs = vapply(tg[name], `[[`, numeric(1), "T_range_obs"),
           T_rep = vapply(tg[name], `[[`, numeric(1), "T_range_rep_mean"),
           absbias = abs(bias))
cat(nrow(s), "fits with a trough check;", sum(s$pass), "passing.\n")
cat("p = P(T_rep >= T_obs) for T_range = peak-to-trough range of OLS-detrended\n",
    "log10(y + 1); p near 1 = the model oscillates too deeply (stage spread too\n",
    "narrow), p near 0 = too shallow, p near 0.5 = calibrated.\n\n", sep = "")

sp <- filter(s, pass)
cat("1. PER CELL. Mean over passing fits of p, of |p - 0.5|, of the observed and\n",
    "replicate-mean peak-to-trough range (log10 units), and mean bias (h) for\n",
    "reference. Unpaired across cells.\n\n", sep = "")
sp |> group_by(window, fit_nc, pin) |>
    summarise(n = n(), p_range = mean(p_range), dev = mean(abs(p_range - 0.5)),
              p_sd = mean(p_sd), T_obs = mean(T_obs), T_rep = mean(T_rep),
              bias = mean(bias), .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)

## Cell-level: is p a monotone read-out of the spread mismatch, and of bias?
vi <- function(b) 1 / (4 * (2 * b + 1))
B <- c(half = 10.0491, truth = 20.0981, x2 = 40.1963)
k168 <- 168 / sp$true_cl[1]
ref_sd <- sqrt(vi(B[["truth"]]) + k168 / 384)
cells <- sp |> group_by(window, fit_nc, pin) |>
    summarise(p = mean(p_range), bias = mean(bias), .groups = "drop") |>
    mutate(spread_sd = if_else(pin == "match", ref_sd,
                               sqrt(vi(unname(B[pin])) + k168 / fit_nc)),
           spread_diff = spread_sd - ref_sd)
cat("\n1b. CELL LEVEL. Spearman over the 11 cells of a window, between the cell's\n",
    "mean p and its spread mismatch at 168 h (cycles, sd) and its mean bias (h);\n",
    "and the spread mismatch at which mean p crosses 0.5 by linear interpolation,\n",
    "which would be 0 for a calibrated check.\n\n", sep = "")
cells |> group_by(window) |>
    summarise(n_cells = n(), rho_p_spread = cor(p, spread_diff, method = "spearman"),
              rho_p_bias = cor(p, bias, method = "spearman"),
              spread_at_p_half = approx(p, spread_diff, xout = 0.5, ties = mean)$y,
              .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 4))) |> as.data.frame() |> print(row.names = FALSE)

cat("\n2. SELECTION. Per (window, rep), candidates = that dataset's passing fits.\n",
    "trough = the candidate with the smallest |p - 0.5|; elpd = highest elpd_loo\n",
    "(_scripts/unitsim-diag.R), where available; oracle = smallest |bias|; random =\n",
    "mean |bias| over the candidates. |bias| in hours.\n\n", sep = "")
dg <- Sys.glob(paste0(D, "diag-*.rds"))
el <- setNames(vapply(dg, \(f) readRDS(f)$elpd_loo, numeric(1)),
               sub("^.*diag-(.*)\\.rds$", "\\1", dg))
sp$elpd <- unname(el[sp$name])
sel <- sp |> group_by(window, rep) |> group_modify(\(g, k) tibble(
    n_cand = nrow(g),
    trough_cell = g$cell[which.min(abs(g$p_range - 0.5))],
    trough_absbias = g$absbias[which.min(abs(g$p_range - 0.5))],
    elpd_absbias = if (all(is.na(g$elpd))) NA_real_ else g$absbias[which.max(g$elpd)],
    ref_absbias = if (any(g$cell == "truth-f384")) g$absbias[g$cell == "truth-f384"] else NA_real_,
    oracle_absbias = min(g$absbias), rand_absbias = mean(g$absbias),
    rho = suppressWarnings(cor(abs(g$p_range - 0.5), g$absbias, method = "spearman")))) |>
    ungroup()
sel |> group_by(window) |>
    summarise(reps = n(), n_cand = mean(n_cand),
              trough = mean(trough_absbias), elpd = mean(elpd_absbias, na.rm = TRUE),
              true_model = mean(ref_absbias, na.rm = TRUE), oracle = mean(oracle_absbias),
              random = mean(rand_absbias),
              trough_beats_elpd = sum(trough_absbias < elpd_absbias, na.rm = TRUE),
              rho = mean(rho, na.rm = TRUE), .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)

cat("\nHow often each cell is chosen by the trough check, per window:\n")
count(sel, window, cell = trough_cell, name = "selected") |> arrange(window, desc(selected)) |>
    as.data.frame() |> print(row.names = FALSE)
