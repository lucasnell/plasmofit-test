# Second simulation study of cycle_length bias: one truth, one fitting
# procedure, data generated at n_c 384 and 1024, always fit at n_c 384. Set up
# 2026-10-21 with Lucas (TODO.md; plan in handoff.md), fixed here BEFORE any
# fit ran. Replaces the 48-cell design of _scripts/simstudy.R, of which 62 of
# 205 fits passed the gate and 9 of 48 cells were readable
# (_data/simstudy-read-2026-10-20-v2.txt): fewer fits, each sampled as a
# production fit.
#
# TRUTH: posterior means of the real-data fit np_wide_total0_nc384_t
#   (wockner-fit.R entry 54; SLURM 30842): no_pool, n_c 384, b_shape
#   estimated per grp_init unit under lognormal(2, 0.5), sd_log10_total0 1,
#   cycle_length prior centred 48 h (sd_logit_cl 1, bounds 35-50 h),
#   adapt_delta 0.95, max_treedepth 12. Mean cycle_length over the 13 trials
#   41.815 h (per trial 39.4-45.6); b_shape 5.86-20.10. The same truth as
#   simstudy.R's truth B at b_shape x 1.
#
# GENERATING n_c, 384 and 1024. The fit is always n_c 384, so gen 384 is
#   correctly specified and gen 1024 is a true n_c 2.7x above the fit's,
#   inside the ML profile's 2-unit interval {512, 1024}
#   (_data/nc-profile-fast-2026-10-09.txt). At gen 1024 every parameter keeps
#   its truth except b_shape, lowered per unit so the TOTAL stage spread at
#   that unit's last observation equals the truth's at 384, since that total
#   is what the data constrain (findings.md, "Why nothing can pin the
#   dispersion"); keeping b_shape would simulate data less dispersed than the
#   real data. Initial age distribution Beta(b, b) over one cycle
#   (beta_starts()), variance 1 / (4 (2b + 1)) cycles^2; the Erlang chain adds
#   k / n_c cycles^2 after k cycles. So
#       1 / (4 (2b' + 1)) = 1 / (4 (2b + 1)) + k_g (1/384 - 1/1024),
#   k_g = unit g's last observation time / its trial's true cycle_length.
#
# FIT, identical to the truth's real-data fit apart from y: production
#   archer_stan_data() defaults as below, n_c 384, b_shape estimated under
#   lognormal(2, 0.5), sd_log10_total0 1; archer_fit(model = "no_pool"),
#   adapt_delta 0.95, max_treedepth 12; 4 chains x (700 warmup + 1000 draws),
#   sampler seed 538065874 (as simstudy.R), chain_id k. No real data is fit.
#
# REPLICATES 1-16 per generating n_c; noise seeds REP_SEEDS below (1-5 are
#   simstudy.R's; 6-16 drawn 2026-10-21 with set.seed(20261021)). The same
#   seed at both gen n_c (common random numbers; the noiseless trajectories
#   differ, so the pairing is partial). Gen-384 data for reps 1-5 are
#   identical to simstudy.R's cell B-g384-f384-bs1.0 (checked in check mode),
#   whose wide arm fitted them with sd_log_b_shape 1.5 and the default sampler.
#   PILOT = reps 1-4 = fits 1-8 = task IDs 1-32 (rep-major order).
#
# LAYOUT: pinned, each chain its own process on 4 physical cores of one L3
#   group (_scripts/simstudy-pinned.sh with SIMSTUDY_RUN / SIMSTUDY_DIR set to
#   this study; CLAUDE.md). 4 threads per chain.
#
# READING RULE, fixed before any fit ran (_scripts/simstudy2-read.R):
#   gate      max R-hat < 1.05 over every monitored quantity (NA excluded),
#             4 chains. A fit that fails is not used for bias.
#   bias      per fit: mean over the 13 trials of (posterior mean of
#             cycle_length[j] - true cycle_length[j]), hours; positive = long.
#   per gen_nc  n passing; mean and se (sd / sqrt(n)) of bias over passing
#             fits; READABLE if >= 3/4 of its replicates pass (12 of 16). Mean
#             over all fits shown as display-only sensitivity, with passing
#             against failing fits.
#   contrast  gen 1024 - gen 384 per rep, both passing (partially paired).
#   b_shape   mean over the 14 units of log(posterior mean / true), against
#             the gen-n_c truth and against the 384 truth.
#   failures  marginal 1.05-1.15, severe >= 1.15; ONE CHAIN APART if a
#             chain's mean lp__ is more than 3 units from the median of the
#             other three.
#   PILOT     reads pass rate, failure bins, and cost per chain. Its bias is
#             printed and labelled display only; it does not decide stage 2.
#
# MODES (SIMSTUDY2_MODE):
#   chain  (default) one chain: task t -> fit ceiling(t / 4), chain
#          (t - 1) %% 4 + 1. Skips if its chain file exists.
#   merge  merge every fit whose 4 chain files exist and that has no summary;
#          write the fit and a one-row summary. Safe to rerun.
#   check  build and simulate fits (SIMSTUDY2_FIT = index or "all"), print
#          the checks listed in the plan; no sampling.
#   SIMSTUDY2_ITER overrides iterations in chain mode for a smoke test, and
#   then writes under a "smoke-" prefix.
#
#   cd /home2/lan68/plasmofit/plasmofit-test
#   SIMSTUDY2_MODE=check SIMSTUDY2_FIT=all Rscript --vanilla _scripts/simstudy2.R
#   sbatch ... _scripts/simstudy-pinned.sh      # see the header there
#   SIMSTUDY2_MODE=merge Rscript --vanilla _scripts/simstudy2.R
#   Rscript --vanilla _scripts/simstudy2-read.R > _data/simstudy2-read-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({
    library(plasmofit)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(rstan)
})

N_CHAINS <- 4L
WARMUP <- 700L
ITER <- 1000L + WARMUP
FIT_NC <- 384L
CONTROL <- list(adapt_delta = 0.95, max_treedepth = 12)
SEED_FIT <- 538065874L                       # as simstudy.R
REP_SEEDS <- c(742160183L, 1908445027L, 355721694L, 1414213562L, 1732050808L,
               273469144L, 923590678L, 832351853L, 1481098644L, 996583259L,
               1100299434L, 938347168L, 681577182L, 1986708724L, 1346550232L,
               1529979755L)
TRUTH <- "np_wide_total0_nc384_t"
THREADS <- as.integer(Sys.getenv("SIMSTUDY2_THREADS", "4"))

## Rep-major, so the pilot (reps 1-4) is fits 1-8 and task IDs 1-32.
FITS <- expand_grid(rep = seq_along(REP_SEEDS), gen_nc = c(384L, 1024L)) |>
    mutate(name = sprintf("g%d-rep%d", gen_nc, rep), fit = row_number())
stopifnot(nrow(FITS) == 32L, !anyDuplicated(FITS$name))

MODE <- Sys.getenv("SIMSTUDY2_MODE", "chain")
iter_env <- Sys.getenv("SIMSTUDY2_ITER", "")
SMOKE <- nzchar(iter_env)
if (SMOKE) { ITER <- as.integer(iter_env); WARMUP <- ITER %/% 2L }
DIR <- "_data/simstudy2/"
PFX <- if (SMOKE) paste0(DIR, "smoke-") else DIR
chain_path <- function(nm, k) sprintf("%schain-%s-c%d.rds", PFX, nm, k)
fit_path <- function(nm) sprintf("%sfit-%s.rds", PFX, nm)
sum_path <- function(nm) sprintf("%ssum-%s.rds", PFX, nm)


# ------------------------------------------------------------------------ #
# Design, truth, and simulated data for one fit
# ------------------------------------------------------------------------ #

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

## Production defaults, as simstudy.R; per-call arguments override them.
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

## The truth's fit was built by wockner-fit.R; its group levels must be the
## ones this script builds, or truth[j] would describe a different unit.
d_truth <- read_rds(sprintf("_data/wock-data-%s.rds", TRUTH))
stopifnot(identical(attr(d_truth, "levels"), attr(d0, "levels")))

read_truth <- function() {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", TRUTH))
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
TRUTH384 <- read_truth()

ends <- cumsum(d0$n_obs); starts <- ends - d0$n_obs + 1L

## Per grp_init unit: its trial (each unit lies in one trial) and the time of
## its last observation in cycles of that trial's true cycle_length.
unit_trial <- vapply(seq_len(d0$n_grp_init), \(g) {
    tr <- unique(d0$grp_cl[d0$grp_init == g]); stopifnot(length(tr) == 1L); tr
}, integer(1))
unit_k_last <- vapply(seq_len(d0$n_grp_init), \(g) {
    s <- which(d0$grp_init == g)
    max(d0$ts[unlist(lapply(s, \(i) starts[i]:ends[i]))]) /
        TRUTH384$cycle_length[unit_trial[g]]
}, numeric(1))

var_init <- function(b) 1 / (4 * (2 * b + 1))          # Beta(b, b), cycles^2
b_from_var <- function(v) (1 / (4 * v) - 1) / 2
total_var <- function(b, n_c) var_init(b) + unit_k_last / n_c

## b_shape at gen_nc that keeps each unit's total spread at its last
## observation equal to the truth's at FIT_NC.
map_b_shape <- function(gen_nc)
    b_from_var(var_init(TRUTH384$b_shape) + unit_k_last * (1 / FIT_NC - 1 / gen_nc))

truth_for <- function(gen_nc) {
    tr <- TRUTH384
    tr$b_shape <- map_b_shape(gen_nc)
    stopifnot(all(tr$b_shape > 2), all(tr$b_shape < d0$max_shape))
    tr
}

## Noiseless trajectories at gen_nc, matrix-exponential path (not the series
## the fit inverts), as simstudy.R.
simulate_mean <- function(truth, gen_nc) {
    y_hat <- numeric(d0$n_total_obs)
    for (s in seq_len(d0$n_ts)) {
        ix <- starts[s]:ends[s]; g <- d0$grp_init[s]
        cl <- truth$cycle_length[d0$grp_cl[s]]
        y0 <- plasmofit:::generate_starts(cl, gen_nc, truth$b_shape[g],
                                          truth$b_offset[g], truth$log10_total0[g])
        y_hat[ix] <- plasmofit:::mat_exp_series(y0, cl, gen_nc,
                                                truth$R[d0$grp_R[s]], d0$mu,
                                                d0$ts[ix], d0$dt_full)
    }
    stopifnot(all(is.finite(y_hat)), all(y_hat >= 0))
    y_hat
}

## Error: normal on log10(y + 1) with sd_iRBC[grp_sd], as the likelihood.
add_noise <- function(y_hat, truth, rep) {
    set.seed(REP_SEEDS[rep])
    sd_obs <- truth$sd_iRBC[rep(d0$grp_sd, d0$n_obs)]
    pmax(0, 10^(log10(y_hat + 1) + rnorm(d0$n_total_obs, 0, sd_obs)) - 1)
}

make_fit_inputs <- function(i) {
    cf <- FITS[i, ]
    truth <- truth_for(cf$gen_nc)
    y_hat <- simulate_mean(truth, cf$gen_nc)
    d <- build(n_c = FIT_NC)
    stopifnot(identical(as.integer(d$grp_cl), as.integer(d0$grp_cl)),
              isTRUE(all.equal(d$ts, d0$ts)), d$n_c == FIT_NC,
              d$b_shape_fixed == 0L)
    d$y <- add_noise(y_hat, truth, cf$rep)
    list(cf = cf, truth = truth, d = d, y_hat = y_hat)
}


# ------------------------------------------------------------------------ #
# Modes
# ------------------------------------------------------------------------ #

dir.create(DIR, showWarnings = FALSE)

if (MODE == "check") {
    which_env <- Sys.getenv("SIMSTUDY2_FIT", "all")
    ids <- if (which_env == "all") FITS$fit else as.integer(which_env)

    ## 1. The fit data are the truth's real-data fit data apart from y (and
    ##    calc_log_lik, which wockner-fit.R set to 1 for loo).
    d <- build(n_c = FIT_NC)
    nm <- setdiff(intersect(names(d), names(d_truth)), c("y", "calc_log_lik"))
    diff_nm <- nm[!vapply(nm, \(n) isTRUE(all.equal(d[[n]], d_truth[[n]])), logical(1))]
    cat("fit data vs the truth's real-data fit data, fields differing (excluding y,",
        "calc_log_lik):", if (length(diff_nm)) paste(diff_nm, collapse = ", ") else "none", "\n")
    stopifnot(length(diff_nm) == 0L)

    ## 2. The b_shape mapping.
    b1024 <- map_b_shape(1024L)
    cat("\nPer grp_init unit: true b_shape at 384 and mapped at 1024; k_last = last",
        "observation in cycles;\nsd_tot = total stage sd (cycles) at k_last;",
        "rel_diff = |sd_tot(1024) / sd_tot(384) - 1|.\n")
    print(data.frame(unit = attr(d0, "levels")[[attr(d0, "grp_columns")$grp_init]],
                     b384 = round(TRUTH384$b_shape, 2), b1024 = round(b1024, 2),
                     k_last = round(unit_k_last, 3),
                     sd_tot384 = round(sqrt(total_var(TRUTH384$b_shape, 384L)), 4),
                     sd_tot1024 = round(sqrt(total_var(b1024, 1024L)), 4),
                     sd_tot1024_unmapped = round(sqrt(total_var(TRUTH384$b_shape, 1024L)), 4),
                     rel_diff = signif(abs(sqrt(total_var(b1024, 1024L) /
                                                total_var(TRUTH384$b_shape, 384L)) - 1), 2)),
          row.names = FALSE)
    stopifnot(all(abs(sqrt(total_var(b1024, 1024L) /
                          total_var(TRUTH384$b_shape, 384L)) - 1) < 1e-10),
              all(b1024 > 2))

    ## 3. Noiseless trajectories: mapped 1024 against 384, and unmapped 1024
    ##    (b_shape kept) against 384, in log10(y + 1). Matched total spread
    ##    should bring the mapped trajectories closer to 384's late in the
    ##    window. Descriptive.
    t0 <- Sys.time(); yh384 <- simulate_mean(truth_for(384L), 384L)
    t384 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    t0 <- Sys.time(); yh1024 <- simulate_mean(truth_for(1024L), 1024L)
    t1024 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    yh1024u <- simulate_mean(TRUTH384, 1024L)
    late <- d0$ts >= 168
    rms <- function(a, b, w) sqrt(mean((log10(a[w] + 1) - log10(b[w] + 1))^2))
    cat(sprintf("\nsimulation time, all 177 series: %.1f s at 384, %.1f s at 1024\n", t384, t1024))
    cat(sprintf("RMS difference from gen 384 in log10(y + 1), noiseless: all obs | obs at >= 168 h\n"))
    cat(sprintf("  gen 1024, b_shape mapped: %.4f | %.4f\n",
                rms(yh1024, yh384, TRUE), rms(yh1024, yh384, late)))
    cat(sprintf("  gen 1024, b_shape kept:   %.4f | %.4f\n",
                rms(yh1024u, yh384, TRUE), rms(yh1024u, yh384, late)))

    ## 4. Each requested fit: simulated data against the real data, and the
    ##    gen-384 reps 1-5 against simstudy.R's data for the same truth.
    cat("\nfit  name        sim mean  sim sd | real mean  real sd | same y as simstudy.R B-g384-f384-bs1.0\n")
    for (i in ids) {
        x <- make_fit_inputs(i)
        same <- "-"
        if (x$cf$gen_nc == 384L && x$cf$rep <= 5L) {
            old <- sprintf("_data/simstudy/chain-B-g384-f384-bs1.0-wide-rep%d-c1.rds", x$cf$rep)
            same <- if (!file.exists(old)) "old chain file missing"
                    else as.character(identical(x$d$y, read_rds(old)$d$y))
        }
        cat(sprintf("%3d  %-10s  %7.4f  %6.4f | %8.4f  %7.4f | %s\n", i, x$cf$name,
                    mean(log10(x$d$y + 1)), sd(log10(x$d$y + 1)),
                    mean(log10(d0$y + 1)), sd(log10(d0$y + 1)), same))
    }
    quit(save = "no")
}

if (MODE == "chain") {
    t <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA))
    stopifnot(!is.na(t), t >= 1, t <= nrow(FITS) * N_CHAINS)
    i <- (t - 1L) %/% N_CHAINS + 1L; k <- (t - 1L) %% N_CHAINS + 1L
    nm <- FITS$name[i]; out <- chain_path(nm, k)
    cat(sprintf("task %d: fit %d (%s), chain %d of %d, %d iter (%d warmup), %d threads, cpus %s\n",
                t, i, nm, k, N_CHAINS, ITER, WARMUP, THREADS,
                sub(".*:\\s*", "", grep("Cpus_allowed_list",
                                        readLines("/proc/self/status"), value = TRUE))))
    if (file.exists(out)) { cat("exists, skipping:", out, "\n"); quit(save = "no") }
    x <- make_fit_inputs(i)
    fk <- archer_fit(x$d, model = "no_pool", chains = 1L, iter = ITER,
                     warmup = WARMUP, seed = SEED_FIT, chain_id = k,
                     threads_per_chain = THREADS, control = CONTROL)
    write_rds(list(fit = fk, d = x$d, truth = x$truth, cf = x$cf, chain = k,
                   seed = SEED_FIT, iter = ITER, warmup = WARMUP,
                   control = CONTROL), out)
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
        cl_mean <- colMeans(as.matrix(f, pars = "cycle_length"))
        bs_mean <- colMeans(as.matrix(f, pars = "b_shape"))
        lp_chain <- apply(rstan::extract(f, "lp__", permuted = FALSE), 2, mean)
        ## largest distance of a chain's mean lp__ from the median of the others
        lp_apart <- max(vapply(seq_along(lp_chain), \(k)
                                   abs(lp_chain[k] - median(lp_chain[-k])), numeric(1)))
        sp <- rstan::get_sampler_params(f, inc_warmup = FALSE)
        lf <- vapply(sp, \(x) sum(x[, "n_leapfrog__"]), numeric(1))
        row <- FITS[i, ] |>
            mutate(max_rhat = max(rh), worst = names(which.max(rh)),
                   n_rhat_over = sum(rh >= 1.05),
                   pass = max(rh) < 1.05,
                   bias = mean(cl_mean - tru$cycle_length),
                   true_cl_mean = mean(tru$cycle_length),
                   est_cl_mean = mean(cl_mean),
                   bs_logratio_gen = mean(log(bs_mean / tru$b_shape)),
                   bs_logratio_384 = mean(log(bs_mean / TRUTH384$b_shape)),
                   lp_spread = diff(range(lp_chain)),
                   lp_apart = lp_apart,
                   div_pct = 100 * mean(unlist(lapply(sp, \(x) x[, "divergent__"]))),
                   treedepth_cap_pct = 100 * mean(unlist(lapply(sp, \(x)
                                           x[, "treedepth__"] >= CONTROL$max_treedepth))),
                   hours_max_chain = max(rowSums(el)) / 3600,
                   hours_sum_chains = sum(el) / 3600)
        row$per_trial <- list(tibble(trial = trial_names, est = cl_mean,
                                     truth = tru$cycle_length))
        row$per_chain <- list(tibble(chain = seq_len(N_CHAINS), lp_mean = lp_chain,
                                     hours = rowSums(el) / 3600,
                                     warmup_hours = el[, "warmup"] / 3600,
                                     leapfrog_per_iter = lf / (ITER - WARMUP),
                                     ms_per_grad = 1e3 * el[, "sample"] / lf))
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

stop("unknown SIMSTUDY2_MODE: ", MODE)
