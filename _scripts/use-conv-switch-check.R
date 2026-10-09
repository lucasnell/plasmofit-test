## Does use_conv change the answer? It must not.
##
## `archer_stan_data(use_conv = 1L)` swaps the Erlang-window polynomial series
## for conv_series() inside the fitted model. The two are the SAME model -- both
## agree with mat_exp_series to ~1e-12 at the function level -- so the posterior
## must not move. This checks that where it matters: the log density and its
## gradient, on the REAL data, at identical parameter values.
##
## This is the gate before any production fit uses use_conv = 1L. A paired
## refit is the stronger test, but it costs hours and should not be started
## until this passes.
##
##
## NOTE, 2026-10-09: conv_series() has been REMOVED from the package (revert
## commit bba580f) because it is slower and less accurate than the
## Erlang-window series at every n_c this project uses. This script therefore
## needs a plasmofit built from commit 05c9c1f to run. It is kept because the
## numbers in claude/findings.md come from it and a claim should name the
## script that produced it; the saved output beside it is the record.
##   cd /home2/lan68/plasmofit/plasmofit-test
##   PLASMOFIT_LIB=/home2/lan68/plasmofit/.Rlib-dev \
##     Rscript --vanilla _scripts/use-conv-switch-check.R \
##     | tee _data/use-conv-switch-check-$(date +%F).txt

LIB <- Sys.getenv("PLASMOFIT_LIB", "/home2/lan68/plasmofit/.Rlib-dev")
if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(rstan); library(readr); library(dplyr)
                                library(plasmofit)})
cat("plasmofit from:", dirname(find.package("plasmofit")), "\n")

## The production grouping, built exactly as _scripts/wockner-fit.R builds it:
## grp_init is trial x inoc_size and grp_sd is trial x cohort. The equality of
## the two forward maps does not depend on the grouping, but using the real one
## means this exercises the same code path a production fit would.
d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = paste(interaction(trial, inoc_size, drop = TRUE)),
           obs_error = paste(interaction(trial, cohort, drop = TRUE)))
cat(sprintf("data: %d observations, %d series, %d trials, %d grp_init, %d grp_sd\n",
            nrow(d), n_distinct(d$id), n_distinct(d$trial),
            n_distinct(d$inoc), n_distinct(d$obs_error)))

mk <- function(uc, n_c) plasmofit::archer_stan_data(
    d, series = "id", time = "time", abundance = "para",
    grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
    grp_sd = "obs_error", n_c = n_c, b_shape = 400, max_shape = 1000,
    sd_log10_total0 = 1, use_conv = uc)

mod <- stan_model(file = file.path(
    dirname(find.package("plasmofit")), "plasmofit", "stan", "archer_fit.stan"),
    model_name = "archer_fit_check")

for (n_c in c(96L, 192L)) {
    sd0 <- mk(0L, n_c); sd1 <- mk(1L, n_c)
    cat(sprintf("\n=== n_c = %d | fft_M would be %d ===\n", n_c,
                plasmofit:::conv_min_M(n_c, max(d$time), sd0$min_cl)))
    f0 <- suppressMessages(sampling(mod, data = sd0, chains = 0))
    f1 <- suppressMessages(sampling(mod, data = sd1, chains = 0))
    set.seed(42)
    U <- lapply(1:5, function(i) rnorm(get_num_upars(f0), 0, 0.5))
    cat("Cells: log posterior density at identical unconstrained parameter\n",
        "values, Erlang-window series against conv_series, on all ", nrow(d),
        "\nobservations. The two are the same model, so this is an equality\n",
        "check. rel is |difference| / |lp|.\n\n", sep = "")
    lp <- t(vapply(U, function(u) {
        a <- log_prob(f0, u); b <- log_prob(f1, u)
        c(ew_series = a, conv_series = b, abs_diff = abs(a - b),
          rel = abs(a - b) / abs(a))
    }, numeric(4)))
    print(as.data.frame(lp), digits = 10)
    g <- t(vapply(U, function(u) {
        a <- grad_log_prob(f0, u); b <- grad_log_prob(f1, u)
        c(max_abs = max(abs(a - b)), rel = max(abs(a - b)) / max(abs(a)))
    }, numeric(2)))
    cat("\ngradient, max over the unconstrained parameters:\n")
    print(as.data.frame(g), digits = 6)
    cat(sprintf("\nworst lp rel %.3e | worst gradient rel %.3e\n",
                max(lp[, "rel"]), max(g[, "rel"])))
    tm <- function(f, u, n = 20) median(replicate(n, system.time(grad_log_prob(f, u))[["elapsed"]]))
    cat(sprintf("gradient cost: ew_series %.4f s | conv_series %.4f s | speedup %.1fx\n",
                tm(f0, U[[1]]), tm(f1, U[[1]]), tm(f0, U[[1]]) / tm(f1, U[[1]])))
    assign(paste0("res", n_c), list(lp = lp, g = g))
}

cat("\n=== verdict ===\n")
worst_lp <- max(res96$lp[, "rel"], res192$lp[, "rel"])
worst_g  <- max(res96$g[, "rel"],  res192$g[, "rel"])
cat(sprintf("worst log_prob relative difference : %.3e\n", worst_lp))
cat(sprintf("worst gradient relative difference : %.3e\n", worst_g))
if (worst_lp < 1e-8 && worst_g < 1e-6) {
    cat("PASS: the switch does not change the model. A paired refit is still the\n",
        "stronger test, but it is now worth the hours.\n", sep = "")
} else {
    cat("FAIL: use_conv changes the posterior. Do not fit with it.\n")
}
