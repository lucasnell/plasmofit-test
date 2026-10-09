## How do the two forward maps scale with n_c, and is there a crossover?
##
## _scripts/use-conv-switch-check.R found conv_series 5x SLOWER than the
## Erlang-window series at n_c = 96, narrowing to 2.8x at 768. The ratio moves
## the right way, so the question is whether it ever reaches 1 at an n_c anyone
## would use. Four points extrapolate to roughly n_c ~ 17,000, but that is an
## extrapolation across a 20-fold gap and the two costs have different shapes:
## the series is ~linear in n_c, the convolution is M log M with M ~ 2 n_c
## rounded up to a power of two, which makes its cost a STAIRCASE rather than a
## smooth curve.
##
## So measure it. The model is already compiled, so changing n_c only changes
## the data and each point costs seconds.
##
## NOTE the large rungs here are computational probes, not biology. n_c is a
## biological assumption in this model and 6144 is nonsense as one; the only
## question being asked is how the two implementations scale.
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
##     Rscript --vanilla _scripts/forward-map-scaling.R \
##     | tee _data/forward-map-scaling-$(date +%F).txt

LIB <- Sys.getenv("PLASMOFIT_LIB", "/home2/lan68/plasmofit/.Rlib-dev")
if (nzchar(LIB)) .libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({library(rstan); library(readr); library(dplyr)
                                library(plasmofit)})
cat("plasmofit from:", dirname(find.package("plasmofit")), "\n")

d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd") |>
    mutate(inoc = paste(interaction(trial, inoc_size, drop = TRUE)),
           obs_error = paste(interaction(trial, cohort, drop = TRUE)))
mk <- function(uc, n_c) plasmofit::archer_stan_data(
    d, series = "id", time = "time", abundance = "para", grp_init = "inoc",
    grp_R = "trial", grp_cl = "trial", grp_sd = "obs_error", n_c = n_c,
    b_shape = 400, max_shape = 1000, sd_log10_total0 = 1, use_conv = uc)
mod <- stan_model(file = file.path(dirname(find.package("plasmofit")),
                                   "plasmofit", "stan", "archer_fit.stan"),
                  model_name = "archer_fit_scaling")

NC <- c(96L, 192L, 384L, 768L, 1536L, 3072L, 6144L)
cat("\nCells: median seconds for one grad_log_prob over all 1130 observations.\n",
    "ew_poly is the Erlang-window series the fitted model uses; conv is\n",
    "conv_series via use_conv = 1. ratio > 1 means the convolution is faster.\n",
    "lp_rel is the relative difference in log_prob, which must stay ~1e-13 or\n",
    "the timing comparison is between two different models.\n\n", sep = "")
res <- NULL
for (n_c in NC) {
    r <- tryCatch({
        f0 <- suppressMessages(sampling(mod, data = mk(0L, n_c), chains = 0))
        f1 <- suppressMessages(sampling(mod, data = mk(1L, n_c), chains = 0))
        set.seed(7); u <- rnorm(get_num_upars(f0), 0, 0.4)
        reps <- if (n_c <= 768) 20 else 5
        tm <- function(f) median(replicate(reps, system.time(grad_log_prob(f, u))[["elapsed"]]))
        a <- tm(f0); b <- tm(f1)
        lp <- c(log_prob(f0, u), log_prob(f1, u))
        data.frame(n_c = n_c, fft_M = plasmofit:::conv_min_M(n_c, max(d$time), 35),
                   ew_poly_s = a, conv_s = b, ratio = a / b,
                   lp_rel = abs(diff(lp)) / abs(lp[1]))
    }, error = function(e) {
        cat(sprintf("n_c = %5d : FAILED -- %s\n", n_c, sub("\n.*", "", conditionMessage(e))))
        NULL
    })
    if (!is.null(r)) { res <- rbind(res, r); print(r, digits = 4, row.names = FALSE) }
    flush.console()
}

cat("\n=== scaling exponents ===\n")
cat("Cells: slope of log(seconds) on log(n_c), i.e. cost ~ n_c^slope, fitted\n",
    "over the rungs measured. A crossover needs the convolution's slope to be\n",
    "clearly smaller.\n\n", sep = "")
fit <- function(y) unname(coef(lm(log(y) ~ log(res$n_c)))[2])
cat(sprintf("  ew_poly   : n_c^%.2f\n", fit(res$ew_poly_s)))
cat(sprintf("  conv      : n_c^%.2f\n", fit(res$conv_s)))
sl <- coef(lm(log(res$ratio) ~ log(res$n_c)))
cat(sprintf("  ratio     : n_c^%.2f\n", unname(sl[2])))
if (unname(sl[2]) > 0) {
    cross <- exp((0 - unname(sl[1])) / unname(sl[2]))
    cat(sprintf("\nextrapolated crossover (ratio = 1) at n_c = %.0f\n", cross))
    cat(if (cross <= 2000)
        "That is inside the range this project would use -- the convolution wins.\n"
        else if (cross <= 10000)
        "That is above any n_c this project would use, but not absurdly so.\n"
        else
        "That is far beyond any usable n_c: the convolution never wins here.\n")
} else {
    cat("\nThe ratio is not improving with n_c; there is no crossover.\n")
}
cat("\nCaveat: an extrapolated crossover is an extrapolation. The measured rungs\n",
    "are what this says anything about.\n", sep = "")
