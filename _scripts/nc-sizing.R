## Thread 5 sizing: what does raising n_c from 96 to 192 cost?
##
## n_c sets the number of sequential exponential compartments, so transit over
## one cycle is Erlang(n_c, n_c / cycle_length) and the stage distribution's sd
## after k cycles is sqrt(k / n_c) cycles. It is currently chosen for numerical
## accuracy by check_erlang_window(), not for biology, which is thread 5.
##
## Before running the real-data arm at n_c = 192 we need its wall time. This
## builds the Wockner data at both values, reports the state dimension and the
## matrix-exponential sizes, then times a SHORT fit at each and reports the
## ratio. It does not produce a usable posterior and makes no scientific claim.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-sizing.R

suppressPackageStartupMessages({
    library(readr); library(dplyr); library(plasmofit)
})

ITER_PROBE  <- as.integer(Sys.getenv("NC_PROBE_ITER",  "20"))
WARMUP_PROBE<- as.integer(Sys.getenv("NC_PROBE_WARMUP","10"))
NC_VALUES   <- c(96L, 192L)

## Same prep as wockner-fit.R. NOTE: no arrange() -- row order sets the group
## factor levels archer_stan_data() derives, and reordering them silently
## invalidated a previous comparison (see claude/gotchas.md).
paras_df <- read_csv(
    if (file.exists("wockner-cleaned.csv")) "wockner-cleaned.csv"
    else "_data/wockner-cleaned.csv",
    col_types = "cccdcdd") |>
    mutate(inoc      = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort,    drop = TRUE) |> paste())

build <- function(nc) {
    archer_stan_data(data = paras_df,
                     series = "id", time = "time", abundance = "para",
                     grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                     grp_sd = "obs_error",
                     n_c = nc,
                     sd_log_b_shape = 1.5,
                     calc_log_lik = 0L)
}

cat("=== build ===\n")
res <- list()
for (nc in NC_VALUES) {
    d <- build(nc)
    cat(sprintf("n_c = %3d | n_ts %d | n_obs %d | dt_full %g | r_dim %s\n",
                d$n_c, d$n_ts, d$n_total_obs, d$dt_full,
                if (!is.null(d$r_dim)) d$r_dim else "NA"))
    res[[as.character(nc)]] <- d
}

cat("\n=== timed probe fit (", ITER_PROBE, " iter, ", WARMUP_PROBE,
    " warmup, 1 chain) ===\n", sep = "")
timing <- list()
for (nc in NC_VALUES) {
    d <- res[[as.character(nc)]]
    t0 <- proc.time()[["elapsed"]]
    f <- archer_fit(d, model = "no_pool", chains = 1L,
                    iter = ITER_PROBE, warmup = WARMUP_PROBE,
                    seed = 1L, threads_per_chain = 1L)
    el <- proc.time()[["elapsed"]] - t0
    lf <- mean(rstan::get_num_leapfrog_per_iteration(f))
    timing[[as.character(nc)]] <- list(elapsed = el, leapfrog = lf)
    cat(sprintf("n_c = %3d | %.1f s wall | %.1f leapfrog/iter | %.4f s per leapfrog\n",
                nc, el, lf, el / (ITER_PROBE * lf)))
}

a <- timing[["96"]]; b <- timing[["192"]]
cat(sprintf("\nper-leapfrog cost ratio 192:96 = %.2fx\n",
            (b$elapsed / (ITER_PROBE * b$leapfrog)) /
            (a$elapsed / (ITER_PROBE * a$leapfrog))))
cat("\nExtrapolate to a production run by multiplying the n_c = 96 arm's known\n",
    "wall time by that ratio. Leapfrog counts differ between a 20-iteration\n",
    "probe and a tuned run, so the PER-LEAPFROG ratio is the transferable\n",
    "number, not the wall-clock ratio.\n", sep = "")
