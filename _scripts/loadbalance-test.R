## How evenly does within-chain threading divide the gradient? (Lucas,
## 2026-10-22: "is there a cheap way to test load balancing too?") No
## sampling: times grad_log_prob only.
##
## The likelihood threads over trajectory combos (reduce_sum over
## traj_combo_partial_sum, one combo per grp_init unit in the full model, 14).
## If combos differ in cost, 4 threads cannot divide the work evenly, and
## threading over (combo, time) pairs instead (a package change) could recover
## the difference. This measures the bound on that.
##
## PART full1 (process pinned to ONE physical core, so one thread):
##   t_c   gradient time with every observation held out except unit c's
##         (hold_out; a combo with no fitted observations skips its trajectory),
##         minus the time with every observation held out (priors only);
##   T1    the full model's gradient time.
## PART full4 (pinned to 4 physical cores of one L3 group, so 4 threads):
##   T4    the full model's gradient time.
## Read together: perfect balance over combos would take max(sum t_c / 4,
## max t_c) + overhead; perfect balance over (combo, time) pairs would take
## sum t_c / 4 + overhead. T4 against those bounds is what finer slicing could
## gain.
## PART unit (one core): within one unit (Mefloquine.1800, 96-216 h, the
##   longest window), the pooled_both gradient time with observations
##   truncated at T = 108, 120, ..., 216 h; increments are the marginal cost of
##   each later time point, i.e. how uneven time points are.
##
## Cost does not depend on y or on parameter values (the Erlang-window bounds
## are fixed from the cycle_length bounds, plasmofit.stan ew_poly_eval), so
## y is the simulated dataset of the single-unit / full-model studies, gen
## n_c 384, rep 1, and the gradient is taken at random unconstrained points.
## Model no_pool (full) / pooled_both (unit), n_c 384, as the studies fit.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   srun --cpus-per-task=8 ... taskset -c <1 core>  Rscript --vanilla _scripts/loadbalance-test.R full1
##   srun --cpus-per-task=8 ... taskset -c <4 cores> Rscript --vanilla _scripts/loadbalance-test.R full4
##   srun ...               taskset -c <1 core>  Rscript --vanilla _scripts/loadbalance-test.R unit
##   (output to _data/loadbalance-<part>-<date>.txt; see _scripts/loadbalance-test.sh)

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(plasmofit); library(readr); library(dplyr); library(rstan)})
part <- commandArgs(TRUE)[1]
stopifnot(part %in% c("full1", "full4", "unit"))
N_EVAL <- 200L
cat("part", part, "| cpus", sub(".*:\\s*", "", grep("Cpus_allowed_list",
    readLines("/proc/self/status"), value = TRUE)), "\n")

paras_df <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste(),
           obs_error = interaction(trial, cohort, drop = TRUE) |> paste())
## Simulated y: the single-unit study's dataset for gen n_c 384, rep 1 is the
## same numbers as simstudy.R's truth-B / gen-384 data (unitsim check); the
## full 1130 observations sit in that study's chain file.
old <- read_rds("_data/simstudy/chain-B-g384-f384-bs1.0-wide-rep1-c1.rds")$d
key_old <- paste(rep(attr(old, "levels")$series, old$n_obs), old$ts)
y_for <- function(d) {
    k <- paste(rep(attr(d, "levels")$series, d$n_obs), d$ts)
    pos <- match(k, key_old); stopifnot(!anyNA(pos)); old$y[pos]
}

## ms per gradient: median of 5 blocks of N_EVAL / 5 evaluations
time_grad <- function(d, model) {
    d$y <- y_for(d)
    mod <- plasmofit:::stanmodels[[model]]
    f <- suppressWarnings(rstan::sampling(mod, data = d, chains = 1, iter = 1, refresh = 0,
                                          seed = 1, algorithm = "Fixed_param"))
    set.seed(11); up <- rnorm(rstan::get_num_upars(f), 0, 0.3)
    invisible(rstan::grad_log_prob(f, up))            # warm-up
    blocks <- vapply(1:5, \(b) {
        t0 <- proc.time()[["elapsed"]]
        for (i in seq_len(N_EVAL / 5)) rstan::grad_log_prob(f, up)
        1e3 * (proc.time()[["elapsed"]] - t0) / (N_EVAL / 5)
    }, numeric(1))
    median(blocks)
}

build_full <- function(df) archer_stan_data(
    data = df, series = "id", time = "time", abundance = "para",
    grp_init = "inoc", grp_R = "trial", grp_cl = "trial", grp_sd = "obs_error",
    calc_log_lik = 0L, n_c = 384L, hold_out = "ho",
    mean_log_b_shape = 2, sd_log_b_shape = 0.5, mean_log10_total0 = 1, sd_log10_total0 = 1)

if (part %in% c("full1", "full4")) {
    units <- unique(paras_df$inoc)
    T_full <- time_grad(build_full(paras_df |> mutate(ho = 0L)), "archer_fit")
    cat(sprintf("full model, all observations fitted: %.3f ms per gradient\n", T_full))
    if (part == "full1") {
        base <- time_grad(build_full(paras_df |> mutate(ho = 1L)), "archer_fit")
        cat(sprintf("all observations held out (priors only): %.3f ms\n\n", base))
        cat("Cells: ms per gradient, one thread, median of 5 blocks of", N_EVAL / 5,
            "evaluations;\nt_c = time with only unit c fitted minus the priors-only time.\n\n")
        res <- do.call(rbind, lapply(units, \(u) {
            df <- paras_df |> mutate(ho = as.integer(inoc != u))
            tu <- time_grad(build_full(df), "archer_fit")
            dfu <- filter(paras_df, inoc == u)
            data.frame(unit = u, series = n_distinct(dfu$id),
                       distinct_times = n_distinct(dfu$time), last_h = max(dfu$time),
                       t_c = round(tu - base, 3))
        }))
        print(res, row.names = FALSE)
        s <- sum(res$t_c)
        cat(sprintf("\nsum t_c %.3f | max t_c %.3f | T1 %.3f | T1 - priors %.3f\n",
                    s, max(res$t_c), T_full, T_full - base))
        cat(sprintf("4-thread bounds: combos max(sum/4, max) = %.3f; (combo, time) pairs sum/4 = %.3f (+ priors %.3f, serial)\n",
                    max(s / 4, max(res$t_c)), s / 4, base))
    }
}

if (part == "unit") {
    u <- "Mefloquine.1800"
    df_u <- filter(paras_df, inoc == u)
    ids <- unique(df_u$id)
    cat("Cells: pooled_both on", u, "alone, observations truncated at T; ms per\n",
        "gradient, one thread, median of 5 blocks; increment = marginal cost of\n",
        "the time points added since the previous T.\n\n")
    cat("(The two-label split computes the unit's trajectory twice per gradient,\n",
        "so on one thread each figure is about twice one trajectory's cost.)\n\n")
    array_fields <- { fd <- build_full(paras_df |> mutate(ho = 0L)); names(fd)[lengths(fd) != 1L] }
    prev <- NA_real_; rows <- list()
    for (Tt in seq(108, 216, by = 12)) {
        dft <- filter(df_u, time <= Tt) |> group_by(id) |> filter(n() >= 2) |> ungroup()
        if (n_distinct(dft$id) < 2) next
        lab <- rep_len(c("A", "B"), n_distinct(dft$id))
        dft$half <- lab[match(dft$id, unique(dft$id))]
        dft$ho <- 0L
        d <- archer_stan_data(data = dft, series = "id", time = "time", abundance = "para",
                              grp_init = "inoc", grp_R = "half", grp_cl = "half",
                              grp_sd = "obs_error", calc_log_lik = 0L, n_c = 384L,
                              hold_out = "ho", mean_log_b_shape = 2, sd_log_b_shape = 0.5,
                              mean_log10_total0 = 1, sd_log10_total0 = 1)
        for (nm in array_fields)
            if (!is.null(d[[nm]]) && length(d[[nm]]) == 1L) d[[nm]] <- as.array(d[[nm]])
        tt <- time_grad(d, "archer_fit_pooled_both")
        rows[[length(rows) + 1]] <- data.frame(T_h = Tt, series = n_distinct(dft$id),
                                               obs = nrow(dft), distinct_times = n_distinct(dft$time),
                                               ms = round(tt, 3), increment = round(tt - prev, 3))
        prev <- tt
    }
    print(do.call(rbind, rows), row.names = FALSE)
}
