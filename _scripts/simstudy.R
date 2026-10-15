# Simulation study of cycle_length bias: how it depends on n_c, on the true
# b_shape, and on how b_shape is fitted. Set up 2026-10-15 (TODO.md,
# "Simulation study of cycle_length bias"); design agreed with Lucas the same
# day, fixed here BEFORE any fit ran.
#
# DESIGN, 48 cells x 5 replicates = 240 fits, 8 chains each:
#
#   truth    A  wock-fit-pooled_cl.rds: posterior means of the n_c 96 pooled_cl
#               fit, one cycle_length (45.012 h) shared by all 13 trials. The
#               truth every earlier simulation used.
#            B  wock-fit-np_wide_total0_nc384_t.rds (entry 54): posterior means
#               of the converged n_c 384 no_pool fit with b_shape estimated,
#               production priors, tuned sampler. PER-TRIAL true cycle lengths
#               (39.4-45.6 h, mean 41.82). Chosen because Lucas ruled out new
#               real-data fits for this study, and it is the only converged
#               n_c 384 fit with production priors and per-unit b_shape (entry
#               43 fixes b_shape at 15 in every unit, so scaling it would give
#               no per-unit variation, and 4 of its 5 runs fail the gate).
#               No converged pooled_cl fit at n_c 384 exists.
#   gen_nc   192, 384        n_c used to GENERATE
#   fit_nc   192, 384        n_c used to FIT (crossed with gen_nc: a 2x2)
#   bs_scale 0.5, 1, 3       true b_shape = the truth fit's per-unit b_shape
#                            times this (b_shape has a lower bound of 2; the
#                            smallest truth is 0.5 x 5.86 = 2.93)
#   method   wide   b_shape estimated, prior lognormal(2, 1.5) on b_shape
#                   (sd_log_b_shape 1.5), the `wide_bshape` prior
#            fixed  b_shape fixed at the simulated truth (diagnostic: only
#                   available in simulation)
#   rep      1-5, noise seeds = wockner-schedule-sim.R's REP_SEEDS[1:5]
#
# Everything else is production: no_pool model, sd_log10_total0 = 1
# (mean 1), cycle_length prior centred at 48 h with sd_logit_cl 1, bounds
# [35, 50], mu = 0, default sampler settings. Simulated data are shared across
# fit_nc and method within (truth, gen_nc, bs_scale, rep), so those contrasts
# are PAIRED on identical data.
#
# Sampling: 8 chains, 700 warmup + 500 post-warmup each (4000 draws, as the
# 4-chain production fits). Each chain is its own SLURM array task pinned to
# ONE logical CPU of one physical core (one Stan thread). Many independent
# fits make throughput, not latency, the target, and the thread benchmark
# (thread-bench-read-2026-10-13.txt) gave the best throughput per core with
# one thread per chain (layout F, 27.1 against 16.0 for one chain per L3
# group).
#
# READING RULE, fixed before any fit ran:
#   gate      max R-hat < 1.05 over every monitored quantity (NA excluded),
#             8 chains. A fit that fails is not used for bias; pass rate per
#             cell is reported as an outcome.
#   bias      per fit: mean over the 13 trials of (posterior mean of
#             cycle_length[j] - true cycle_length[j]), hours; positive = long.
#   per cell  mean, and range, of bias over the replicates that pass;
#             n passing of 5. A cell with fewer than 3 passing replicates is
#             reported as NOT READABLE, not as a bias.
#   contrasts paired on the same simulated data, per replicate, both fits
#             passing: wide - fixed (cost of estimating b_shape), and fit_nc
#             384 - 192 at each gen_nc (misspecification = off-diagonal minus
#             diagonal).
#   also      b_shape recovery in the wide arm: mean over the 14 units of
#             log(posterior mean b_shape / true b_shape).
#
# MODES (SIMSTUDY_MODE):
#   chain  (default) one chain: task t -> fit ceiling(t / 8), chain
#          (t - 1) %% 8 + 1. Skips if its chain file exists.
#   merge  merge every fit whose 8 chain files exist and that has no merged
#          file; write the fit and a one-row summary. Safe to rerun.
#   check  build and simulate one fit (SIMSTUDY_FIT = index), print, no
#          sampling. SIMSTUDY_ITER overrides iterations in chain mode for a
#          smoke test, and then writes under a "smoke-" prefix.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   sbatch _scripts/simstudy.sh            # 1920 chain tasks
#   SIMSTUDY_MODE=merge Rscript --vanilla _scripts/simstudy.R
#   Rscript --vanilla _scripts/simstudy-read.R > _data/simstudy-read-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({
    library(plasmofit)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(rstan)
})

N_CHAINS <- 8L
WARMUP <- 700L
ITER <- 500L + WARMUP
SEED_FIT <- 538065874L                       # as wockner-schedule-sim.R
REP_SEEDS <- c(742160183L, 1908445027L, 355721694L, 1414213562L, 1732050808L)
TRUTHS <- c(A = "pooled_cl", B = "np_wide_total0_nc384_t")

## Rep-major, so replicate 1 of every cell finishes first; within a rep the
## n_c 384 fits first, since they are the longest.
FITS <- expand_grid(rep = seq_along(REP_SEEDS),
                    fit_nc = c(384L, 192L),
                    truth = names(TRUTHS),
                    gen_nc = c(192L, 384L),
                    bs_scale = c(0.5, 1, 3),
                    method = c("wide", "fixed")) |>
    mutate(cell = sprintf("%s-g%d-f%d-bs%s-%s", truth, gen_nc, fit_nc,
                          format(bs_scale), method),
           name = sprintf("%s-rep%d", cell, rep),
           fit = row_number())
stopifnot(nrow(FITS) == 240L, !anyDuplicated(FITS$name))

MODE <- Sys.getenv("SIMSTUDY_MODE", "chain")
iter_env <- Sys.getenv("SIMSTUDY_ITER", "")
SMOKE <- nzchar(iter_env)
if (SMOKE) { ITER <- as.integer(iter_env); WARMUP <- ITER %/% 2L }
PFX <- if (SMOKE) "_data/simstudy/smoke-" else "_data/simstudy/"
chain_path <- function(nm, k) sprintf("%schain-%s-c%d.rds", PFX, nm, k)
fit_path <- function(nm) sprintf("%sfit-%s.rds", PFX, nm)
sum_path <- function(nm) sprintf("%ssum-%s.rds", PFX, nm)


# ------------------------------------------------------------------------ #
# Design, truth, and simulated data for one fit
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

## Production defaults; per-arm arguments override them.
build <- function(...) {
    args <- modifyList(list(data = paras_df, series = "id", time = "time",
                            abundance = "para", grp_init = "inoc",
                            grp_R = "trial", grp_cl = "trial",
                            grp_sd = "obs_error", calc_log_lik = 0L,
                            mean_log_b_shape = 2, sd_log_b_shape = 0.5,
                            mean_log10_total0 = 1, sd_log10_total0 = 1),
                       list(...))
    do.call(archer_stan_data, args)
}
d0 <- build()

## Truth B's fit was built by wockner-fit.R; its group levels must be the
## ones this script builds, or truth[j] would describe a different unit.
lv_B <- attr(read_rds("_data/wock-data-np_wide_total0_nc384_t.rds"), "levels")
stopifnot(identical(lv_B, attr(d0, "levels")))

read_truth <- function(tr) {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", TRUTHS[[tr]]))
    pm <- function(p) unname(colMeans(as.matrix(f, pars = p)))
    out <- list(b_shape = pm("b_shape"), b_offset = pm("b_offset"),
                log10_total0 = pm("log10_total0"), R = pm("R"),
                sd_iRBC = pm("sd_iRBC"), cycle_length = pm("cycle_length"))
    stopifnot(length(out$b_shape) == d0$n_grp_init,
              length(out$b_offset) == d0$n_grp_init,
              length(out$log10_total0) == d0$n_grp_init,
              length(out$R) == d0$n_grp_R, length(out$sd_iRBC) == d0$n_grp_sd,
              length(out$cycle_length) == d0$n_grp_cl,
              all(out$cycle_length > d0$min_cl), all(out$cycle_length < d0$max_cl))
    out
}

make_fit_inputs <- function(i) {
    cf <- FITS[i, ]
    truth <- read_truth(cf$truth)
    truth$b_shape <- truth$b_shape * cf$bs_scale
    stopifnot(all(truth$b_shape > 2), all(truth$b_shape < d0$max_shape))

    ## Simulate at gen_nc with the matrix-exponential path, not the series the
    ## fit inverts, as wockner-schedule-sim.R does. Error: normal on
    ## log10(y + 1) with sd_iRBC[grp_sd], as the likelihood.
    ends <- cumsum(d0$n_obs); starts <- ends - d0$n_obs + 1L
    y_hat <- numeric(d0$n_total_obs)
    for (s in seq_len(d0$n_ts)) {
        ix <- starts[s]:ends[s]; g <- d0$grp_init[s]
        cl <- truth$cycle_length[d0$grp_cl[s]]
        y0 <- plasmofit:::generate_starts(cl, cf$gen_nc, truth$b_shape[g],
                                          truth$b_offset[g], truth$log10_total0[g])
        y_hat[ix] <- plasmofit:::mat_exp_series(y0, cl, cf$gen_nc,
                                                truth$R[d0$grp_R[s]], d0$mu,
                                                d0$ts[ix], d0$dt_full)
    }
    stopifnot(all(is.finite(y_hat)), all(y_hat >= 0))
    set.seed(REP_SEEDS[cf$rep])
    sd_obs <- truth$sd_iRBC[rep(d0$grp_sd, d0$n_obs)]
    y_sim <- pmax(0, 10^(log10(y_hat + 1) + rnorm(d0$n_total_obs, 0, sd_obs)) - 1)

    d <- if (cf$method == "fixed") build(n_c = cf$fit_nc, b_shape = truth$b_shape)
         else build(n_c = cf$fit_nc, sd_log_b_shape = 1.5)
    stopifnot(identical(as.integer(d$grp_cl), as.integer(d0$grp_cl)),
              isTRUE(all.equal(d$ts, d0$ts)), d$n_c == cf$fit_nc)
    d$y <- y_sim
    list(cf = cf, truth = truth, d = d, y_hat = y_hat)
}


# ------------------------------------------------------------------------ #
# Modes
# ------------------------------------------------------------------------ #

dir.create("_data/simstudy", showWarnings = FALSE)

if (MODE == "check") {
    i <- as.integer(Sys.getenv("SIMSTUDY_FIT", "1"))
    x <- make_fit_inputs(i)
    print(as.data.frame(x$cf))
    cat(sprintf("true cycle_length %.3f-%.3f h (mean %.3f); true b_shape %.2f-%.2f\n",
                min(x$truth$cycle_length), max(x$truth$cycle_length),
                mean(x$truth$cycle_length), min(x$truth$b_shape), max(x$truth$b_shape)))
    cat(sprintf("noiseless vs real log10(y+1): r = %.3f; sim mean %.3f sd %.3f; real mean %.3f sd %.3f\n",
                cor(log10(x$y_hat + 1), log10(d0$y + 1)),
                mean(log10(x$d$y + 1)), sd(log10(x$d$y + 1)),
                mean(log10(d0$y + 1)), sd(log10(d0$y + 1))))
    cat("fit n_c", x$d$n_c, "| b_shape fixed:", x$d$b_shape_fixed,
        "| sd_log_b_shape", unique(x$d$sd_log_b_shape),
        "| sd_log10_total0", unique(x$d$sd_log10_total0), "\n")
    quit(save = "no")
}

if (MODE == "chain") {
    t <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
    stopifnot(!is.na(t), t >= 1, t <= nrow(FITS) * N_CHAINS)
    i <- (t - 1L) %/% N_CHAINS + 1L; k <- (t - 1L) %% N_CHAINS + 1L
    nm <- FITS$name[i]; out <- chain_path(nm, k)
    cat(sprintf("task %d: fit %d (%s), chain %d of %d, %d iter (%d warmup), cpus %s\n",
                t, i, nm, k, N_CHAINS, ITER, WARMUP,
                sub(".*:\\s*", "", grep("Cpus_allowed_list",
                                        readLines("/proc/self/status"), value = TRUE))))
    if (file.exists(out)) { cat("exists, skipping:", out, "\n"); quit(save = "no") }
    x <- make_fit_inputs(i)
    fk <- archer_fit(x$d, model = "no_pool", chains = 1L, iter = ITER,
                     warmup = WARMUP, seed = SEED_FIT, chain_id = k,
                     threads_per_chain = 1L)
    write_rds(list(fit = fk, d = x$d, truth = x$truth, cf = x$cf, chain = k,
                   seed = SEED_FIT, iter = ITER, warmup = WARMUP), out)
    print(rstan::get_elapsed_time(fk))
    cat("wrote", out, "\n")
    quit(save = "no")
}

if (MODE == "merge") {
    trial_names <- attr(d0, "levels")[[attr(d0, "grp_columns")$grp_cl]]
    n_done <- 0L
    for (i in seq_len(nrow(FITS))) {
        nm <- FITS$name[i]
        paths <- vapply(seq_len(N_CHAINS), \(k) chain_path(nm, k), character(1))
        if (file.exists(sum_path(nm)) || !all(file.exists(paths))) next
        parts <- lapply(paths, read_rds)
        ## Every chain must have fitted the same data under the same truth.
        for (pt in parts) stopifnot(identical(pt$d, parts[[1]]$d),
                                    identical(pt$truth, parts[[1]]$truth),
                                    pt$seed == SEED_FIT)
        stopifnot(identical(vapply(parts, `[[`, integer(1), "chain"), seq_len(N_CHAINS)))
        f <- rstan::sflist2stanfit(lapply(parts, `[[`, "fit"))
        tru <- parts[[1]]$truth
        el <- do.call(rbind, lapply(parts, \(p) rstan::get_elapsed_time(p$fit)))
        rm(parts)
        write_rds(f, fit_path(nm))

        s <- rstan::summary(f)$summary
        rh <- s[, "Rhat"]; rh <- rh[is.finite(rh)]
        cl <- as.matrix(f, pars = "cycle_length")
        cl_mean <- colMeans(cl)
        bs_mean <- if (FITS$method[i] == "wide") colMeans(as.matrix(f, pars = "b_shape")) else NA
        lp <- rstan::extract(f, "lp__", permuted = FALSE)
        sp <- rstan::get_sampler_params(f, inc_warmup = FALSE)
        row <- FITS[i, ] |>
            mutate(max_rhat = max(rh), worst = names(which.max(rh)),
                   n_rhat_over = sum(rh >= 1.05),
                   pass = max(rh) < 1.05,
                   bias = mean(cl_mean - tru$cycle_length),
                   true_cl_mean = mean(tru$cycle_length),
                   est_cl_mean = mean(cl_mean),
                   bs_logratio = if (all(is.na(bs_mean))) NA_real_
                                 else mean(log(bs_mean / tru$b_shape)),
                   div_pct = 100 * mean(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                   treedepth_cap_pct = 100 * mean(unlist(lapply(sp, \(x)
                                           x[, "treedepth__"] >= 10))),
                   lp_spread = diff(range(apply(lp, 2, mean))),
                   hours_max_chain = max(rowSums(el)) / 3600,
                   hours_sum_chains = sum(el) / 3600)
        row$per_trial <- list(tibble(trial = trial_names, est = cl_mean,
                                     truth = tru$cycle_length))
        write_rds(row, sum_path(nm))
        n_done <- n_done + 1L
        cat(sprintf("merged %s: max R-hat %.3f %s, bias %+.3f h, %.2f h longest chain\n",
                    nm, row$max_rhat, if (row$pass) "PASS" else "FAIL",
                    row$bias, row$hours_max_chain))
    }
    cat("merged", n_done, "fits this run;",
        length(Sys.glob(paste0(PFX, "sum-*.rds"))), "of", nrow(FITS), "now merged\n")
    quit(save = "no")
}

stop("unknown SIMSTUDY_MODE: ", MODE)
