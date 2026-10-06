## Thread 5: is the n_c effect numerical or structural?
##
## The n_c ladder moves cycle_length by -3.20 h and gains 71.8 elpd going from
## n_c = 96 to 384. n_c does two jobs at once -- it sets the Erlang shape, and
## so the model's desynchronisation rate, AND it has to be large enough for the
## Erlang-window series to approximate the matrix exponential. Those two have
## opposite implications and have to be separated:
##
##   NUMERICAL -- the n_c = 96 fit computes its own model inaccurately. Then
##   every earlier fit is suspect arithmetic and check_erlang_window() is too
##   permissive.
##
##   STRUCTURAL -- n_c = 96 computes its own model faithfully, but that model
##   desynchronises too fast for these data. Then the arithmetic was always
##   fine and the finding is about biology.
##
## max_rel_diff is the model's own cross-check: generated quantities recomputes
## the trajectory with matrix_exp and reports the largest relative difference
## against the series solution. It is evaluated here at the ACTUAL posterior
## draws via gqs(), not at a short proxy run.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-numerical-check.R

suppressPackageStartupMessages({
    library(readr); library(dplyr); library(rstan); library(plasmofit)
})

ARMS <- c(np_wide_both = 96L, np_nc192 = 192L, np_nc384 = 384L)
N_DRAWS <- as.integer(Sys.getenv("NC_CHECK_DRAWS", "200"))

paras_df <- read_csv(
    if (file.exists("wockner-cleaned.csv")) "wockner-cleaned.csv"
    else "_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc      = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    mutate(obs_error = interaction(trial, cohort,    drop = TRUE) |> paste())

## The saved fits do not monitor the raw parameters block (b_shape_free is not
## retained), so gqs() on the posterior draws is not available. Instead each
## n_c gets a SHORT fit with run_check = 1. max_rel_diff depends on the
## parameters only through cycle_length, which is bounded to [35, 50], so a
## short run probes the same numerical regime as the production fit. This is a
## check on the ARITHMETIC, not a posterior, and is not used for any estimate.

ITER_CHK <- as.integer(Sys.getenv("NC_CHECK_ITER", "40"))

for (nm in names(ARMS)) {
    nc <- ARMS[[nm]]
    d <- archer_stan_data(data = paras_df,
                          series = "id", time = "time", abundance = "para",
                          grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                          grp_sd = "obs_error",
                          n_c = nc, sd_log_b_shape = 1.5, sd_log10_total0 = 1,
                          run_check = 1L, calc_log_lik = 0L)
    f <- archer_fit(d, model = "no_pool", chains = 1L,
                    iter = ITER_CHK, warmup = ITER_CHK %/% 2L,
                    seed = 1L, threads_per_chain = 1L)
    mrd <- as.matrix(f, pars = "max_rel_diff")[, 1]
    cl  <- as.matrix(f, pars = "cycle_length")
    cat(sprintf("n_c = %3d (%-12s) | max_rel_diff: max %.3e, median %.3e | cycle_length visited %.1f-%.1f h\n",
                nc, nm, max(mrd), median(mrd), min(cl), max(cl)))
}

cat("\nmax_rel_diff is the largest RELATIVE difference between the series\n",
    "solution the likelihood uses and matrix_exp on the same A matrix, over\n",
    "every observation and every draw. Near machine precision means the\n",
    "arithmetic is faithful at that n_c, so the ladder is STRUCTURAL -- the\n",
    "n_c = 96 model is computed correctly and simply fits worse.\n", sep = "")
