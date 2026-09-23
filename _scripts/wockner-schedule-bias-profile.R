# Where does the +1.5 h cycle_length bias come from, and how does it scale?
#
# wockner-schedule-sim.R found that fitting no_pool to data simulated from a
# single known cycle_length recovers it ~1.5 h too high, with the prior's
# influence removed (see claude/CLAUDE.md, "Schedule-bias simulation"). That
# result came from full hierarchical fits, so prior, hierarchy, shrinkage and
# likelihood are all still tangled together in it.
#
# This strips everything away but the likelihood. For each trial's real
# observation schedule, simulate from the known truth and find the maximum
# likelihood cycle_length directly, jointly with that trial's own nuisance
# parameters. No prior, no hierarchy, no MCMC: if the maximum sits above the
# truth, the observation schedules bias the likelihood itself, which is the
# one explanation no change of prior can fix.
#
# It is cheap because all series within a trial share grp_cl and grp_R, and
# almost always one grp_init, so they share a single trajectory and differ
# only in which times they observe -- the same dedup archer_fit.stan does. One
# full-grid solve per likelihood evaluation (~15 ms), not one per series.
#
# The control matters as much as the result: with noiseless data the maximum
# must land on the truth. If it does not, the finding is a bug in this script
# rather than a property of the design, and nothing below can be read.
#
#   srun -N 1 -n 1 -c 8 --mem=16G --time=02:00:00 \
#       Rscript --vanilla _scripts/wockner-schedule-bias-profile.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(parallel)
})

N_REP <- 16L                  # noise realizations per trial
SEED <- 605114982L
N_CORE <- max(1L, min(8L, detectCores()))


# ------------------------------------------------------------------------ #
# The real design and the same truth the simulation used
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- archer_stan_data(paras_df,
                      series = "id", time = "time", abundance = "para",
                      grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                      grp_sd = "obs_error", calc_log_lik = 0L)
trial_names <- attr(d, "levels")[[attr(d, "grp_columns")$grp_cl]]

post_mean <- function(f, par) unname(colMeans(as.matrix(f, pars = par)))
f_pl <- read_rds("_data/wock-fit-pooled_cl.rds")
truth <- list(b_shape      = post_mean(f_pl, "b_shape"),
              b_offset     = post_mean(f_pl, "b_offset"),
              log10_total0 = post_mean(f_pl, "log10_total0"),
              R            = post_mean(f_pl, "R"),
              sd_iRBC      = post_mean(f_pl, "sd_iRBC"))
truth$cycle_length <- post_mean(f_pl, "cycle_length")[1]
rm(f_pl); invisible(gc())

ends <- cumsum(d$n_obs); starts <- ends - d$n_obs + 1L
n_full <- as.integer(round(max(d$ts) / d$dt_full)) + 1L   ## matches archer_fit.stan:179

## everything one trial needs, gathered once
trial_data <- map(seq_len(d$n_grp_cl), \(j) {
    sers <- which(d$grp_cl == j)
    list(j = j,
         name = trial_names[j],
         gi = d$grp_init[sers],
         gi_set = sort(unique(d$grp_init[sers])),
         idx = map(sers, \(i) as.integer(d$ts[starts[i]:ends[i]] / d$dt_full) + 1L),
         sd = truth$sd_iRBC[d$grp_sd[sers]],
         n_obs = d$n_obs[sers],
         R = truth$R[d$grp_R[sers[1]]])
})


# ------------------------------------------------------------------------ #
# Likelihood
# ------------------------------------------------------------------------ #

lgt <- function(p) log(p / (1 - p))
ilgt <- function(x) 1 / (1 + exp(-x))
cl_of <- function(z) d$min_cl + (d$max_cl - d$min_cl) * ilgt(z)

## par = c(z_cl, log_b_shape[G], z_b_offset[G], log10_total0[G], log_R)
unpack <- function(par, G) {
    list(cl = cl_of(par[1]),
         b_shape = exp(par[1 + seq_len(G)]),
         b_offset = ilgt(par[1 + G + seq_len(G)]),
         log10_total0 = par[1 + 2 * G + seq_len(G)],
         R = exp(par[2 + 3 * G]))
}

y_hat_trial <- function(p, td) {
    G <- length(td$gi_set)
    traj <- vector("list", G)
    for (g in seq_len(G)) {
        y0 <- plasmofit:::generate_starts(p$cl, d$n_c, p$b_shape[g],
                                          p$b_offset[g], p$log10_total0[g])
        traj[[g]] <- plasmofit:::full_mat_exp_series(y0, p$cl, d$n_c, p$R, d$mu,
                                                     d$dt_full, n_full)
    }
    gmap <- match(td$gi, td$gi_set)
    map2(td$idx, gmap, \(ix, g) traj[[g]][ix])
}

nll <- function(par, td, y_obs_trans) {
    G <- length(td$gi_set)
    p <- unpack(par, G)
    if (!all(is.finite(unlist(p))) || p$b_shape[1] > d$max_shape ||
        p$R > d$max_R) return(1e12)
    yh <- tryCatch(y_hat_trial(p, td), error = function(e) NULL)
    if (is.null(yh)) return(1e12)
    mu_trans <- log10(unlist(yh) + 1)
    if (!all(is.finite(mu_trans))) return(1e12)
    sd_v <- rep(td$sd, td$n_obs)
    -sum(dnorm(y_obs_trans, mu_trans, sd_v, log = TRUE))
}

par_at <- function(td, cl) {
    g <- td$gi_set
    c(lgt((cl - d$min_cl) / (d$max_cl - d$min_cl)),
      log(truth$b_shape[g]), lgt(truth$b_offset[g]),
      truth$log10_total0[g], log(td$R))
}

## Maximize from several cycle_length starts and keep the best, so an aliased
## secondary mode is found rather than silently settled into.
fit_mle <- function(td, y_obs_trans, cl_starts) {
    best <- NULL
    for (cl0 in cl_starts) {
        o <- tryCatch(optim(par_at(td, cl0), nll, td = td,
                            y_obs_trans = y_obs_trans, method = "Nelder-Mead",
                            control = list(maxit = 4000, reltol = 1e-11)),
                      error = function(e) NULL)
        if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    if (is.null(best)) return(NA_real_)
    cl_of(best$par[1])
}

CL_STARTS <- truth$cycle_length + c(-2.5, 0, 2.5)


# ------------------------------------------------------------------------ #
# Control: noiseless data. The maximum must land on the truth.
# ------------------------------------------------------------------------ #

cat("=== control: noiseless data, maximum should equal the truth ===\n")
cat(sprintf("  true cycle_length %.4f h\n", truth$cycle_length))

ctrl <- mclapply(trial_data, \(td) {
    p <- unpack(par_at(td, truth$cycle_length), length(td$gi_set))
    yt <- log10(unlist(y_hat_trial(p, td)) + 1)
    fit_mle(td, yt, CL_STARTS)
}, mc.cores = N_CORE) |> unlist()

ctrl_tab <- tibble(trial = trial_names, mle = ctrl,
                   err = ctrl - truth$cycle_length)
print(ctrl_tab |> mutate(across(where(is.numeric), \(x) round(x, 5))), n = Inf)
cat(sprintf("\n  max |error|: %.2e h\n", max(abs(ctrl_tab$err))))
if (max(abs(ctrl_tab$err)) > 0.01) {
    cat("\n  CONTROL FAILED. The optimizer does not recover the truth from\n",
        "  noiseless data, so nothing below is interpretable.\n", sep = "")
} else {
    cat("  Control passes: the estimator is unbiased without noise, so any\n",
        "  bias below is a noise-by-schedule effect.\n", sep = "")
}


# ------------------------------------------------------------------------ #
# With noise, at each trial's real schedule
# ------------------------------------------------------------------------ #

set.seed(SEED)
seeds <- sample.int(.Machine$integer.max, N_REP)

res <- mclapply(seq_len(N_REP), \(r) {
    set.seed(seeds[r])
    ## draw every trial's noise before optimizing, so replicate r is a
    ## coherent dataset rather than 13 independent ones
    yobs <- map(trial_data, \(td) {
        p <- unpack(par_at(td, truth$cycle_length), length(td$gi_set))
        yh <- unlist(y_hat_trial(p, td))
        sd_v <- rep(td$sd, td$n_obs)
        log10(pmax(0, 10^(log10(yh + 1) + rnorm(length(yh), 0, sd_v)) - 1) + 1)
    })
    map2_dbl(trial_data, yobs, \(td, yt) fit_mle(td, yt, CL_STARTS))
}, mc.cores = N_CORE)

mle <- do.call(rbind, res)            # N_REP x n_trial
colnames(mle) <- trial_names

meta <- paras_df |>
    summarise(.by = trial, n_series = n_distinct(id), n_obs = n()) |>
    mutate(obs_per_series = n_obs / n_series)
ops <- meta$obs_per_series[match(trial_names, meta$trial)]

out <- tibble(trial = trial_names,
              obs_per_series = ops,
              mle_mean = colMeans(mle, na.rm = TRUE),
              mle_sd = apply(mle, 2, sd, na.rm = TRUE),
              bias = colMeans(mle, na.rm = TRUE) - truth$cycle_length) |>
    arrange(obs_per_series)

cat("\n=== maximum-likelihood cycle_length, ", N_REP,
    " noise replicates per trial ===\n", sep = "")
out |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf)

cat(sprintf("\n  mean bias across trials: %+.3f h\n", mean(out$bias)))
cat(sprintf("  corr(bias, obs_per_series): %+.3f\n",
            cor(out$bias, out$obs_per_series)))

## The same correlation the full fits report, but computed on pure MLEs: no
## prior, no hierarchy, no shrinkage.
r_rep <- apply(mle, 1, \(x) cor(x, ops))
cat(sprintf("  corr(MLE, obs_per_series) within a replicate: mean %+.3f (95%% %+.3f to %+.3f)\n",
            mean(r_rep), quantile(r_rep, 0.025), quantile(r_rep, 0.975)))
cat("  real data, posterior means: -0.933\n")



# ------------------------------------------------------------------------ #
# The pooled likelihood: one cycle_length shared by all 13 trials
#
# The hierarchical fits' bias is in a POPULATION estimate, which pools 13
# trials, so the average of 13 per-trial MLEs above is not the matching
# quantity. This is: a single cycle_length shared by every trial, with each
# trial keeping its own nuisance parameters. No prior, no hierarchy.
#
# Given cycle_length the trials are separable, so this profiles over a grid of
# shared cycle_length and optimizes each trial's nuisance parameters inside,
# warm-started along the grid. If the pooled maximum lands on the truth, the
# bias in the hierarchical fits comes from the hierarchy and the prior, not
# from the likelihood at these schedules.
# ------------------------------------------------------------------------ #

N_REP_POOL <- 8L
CL_GRID <- seq(41, 49, by = 0.5)

prof_rep <- function(r) {
    set.seed(seeds[r])
    yobs <- map(trial_data, \(td) {
        p <- unpack(par_at(td, truth$cycle_length), length(td$gi_set))
        yh <- unlist(y_hat_trial(p, td))
        sd_v <- rep(td$sd, td$n_obs)
        log10(pmax(0, 10^(log10(yh + 1) + rnorm(length(yh), 0, sd_v)) - 1) + 1)
    })
    ## nuisance-only objective at a fixed cycle_length
    nll_fix <- function(q, td, yt, cl) {
        z <- lgt((cl - d$min_cl) / (d$max_cl - d$min_cl))
        nll(c(z, q), td, yt)
    }
    warm <- map(trial_data, \(td) par_at(td, truth$cycle_length)[-1])
    tot <- numeric(length(CL_GRID))
    for (k in seq_along(CL_GRID)) {
        cl <- CL_GRID[k]
        v <- 0
        for (j in seq_along(trial_data)) {
            o <- tryCatch(optim(warm[[j]], nll_fix, td = trial_data[[j]],
                                yt = yobs[[j]], cl = cl, method = "Nelder-Mead",
                                control = list(maxit = 3000, reltol = 1e-10)),
                          error = function(e) NULL)
            if (!is.null(o)) { warm[[j]] <- o$par; v <- v + o$value }
            else v <- v + 1e12
        }
        tot[k] <- v
    }
    tot
}

prof <- mclapply(seq_len(N_REP_POOL), prof_rep, mc.cores = N_CORE)
prof_mat <- do.call(rbind, prof)

## parabolic refinement around the grid minimum
argmin_refine <- function(v) {
    k <- which.min(v)
    if (k == 1L || k == length(v)) return(CL_GRID[k])
    y1 <- v[k - 1]; y2 <- v[k]; y3 <- v[k + 1]
    h <- CL_GRID[2] - CL_GRID[1]
    CL_GRID[k] + h * (y1 - y3) / (2 * (y1 - 2 * y2 + y3))
}
pooled_mle <- apply(prof_mat, 1, argmin_refine)

cat("\n=== pooled likelihood: one cycle_length shared by all 13 trials ===\n")
cat(sprintf("  %d replicates, grid %.0f-%.0f by %.1f h\n",
            N_REP_POOL, min(CL_GRID), max(CL_GRID), CL_GRID[2] - CL_GRID[1]))
cat(sprintf("  pooled MLE: %s\n", paste(sprintf("%.2f", pooled_mle), collapse = " ")))
cat(sprintf("  mean %.3f h, sd %.3f, bias %+.3f h (truth %.3f)\n",
            mean(pooled_mle), sd(pooled_mle),
            mean(pooled_mle) - truth$cycle_length, truth$cycle_length))
cat(sprintf("  se of the mean bias: %.3f h\n", sd(pooled_mle) / sqrt(N_REP_POOL)))
cat("\n  Compare against the hierarchical fits' population estimates, which\n",
    "  ran +1.0 to +2.6 h high on the same design.\n", sep = "")

write_rds(list(truth = truth$cycle_length, control = ctrl_tab, mle = mle,
               obs_per_series = ops, summary = out, r_rep = r_rep,
               cl_grid = CL_GRID, profile = prof_mat, pooled_mle = pooled_mle),
          "_data/wock-schedbias-profile.rds")
cat("\nwrote _data/wock-schedbias-profile.rds\n")
