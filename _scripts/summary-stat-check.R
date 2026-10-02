# Would reporting the posterior MEDIAN instead of the mean reduce the
# cycle-length bias, and where else would it change a number?
#
# Two different questions get conflated here and only one of them is about
# skew:
#
#   1. SUMMARY CHOICE. mean vs median vs mode of a marginal posterior. This
#      matters exactly as much as that marginal is skewed, and no more.
#   2. MARGINALISATION. The posterior mean, median and marginal mode are ALL
#      marginal summaries -- they integrate over the other ~200 parameters
#      alike. If the bias comes from integrating over a ridge and two
#      unidentified nuisances, every one of them carries it equally and
#      switching between them cannot help. Only a JOINT summary (the MAP)
#      escapes that, which is why it is the discriminating test.
#
# This script answers 1 with numbers. It cannot answer 2.
#
# Simulation fits have a known truth, so bias is measurable; the real fit is
# included to show where the summary choice changes a REPORTED number even
# though no bias can be computed for it.
#
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/summary-stat-check.R 2>&1 | tee _data/summary-stat-check.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

TRUE_CL <- 45.012   # the value the `default` arm's data were simulated from

skew <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3

## per-group mean / median / skew for one parameter of one fit
grab <- function(f, par) {
    m <- as.matrix(f, pars = par)
    tibble(par = par, idx = seq_len(ncol(m)),
           mean = colMeans(m),
           median = apply(m, 2, stats::median),
           skew = apply(m, 2, skew))
}

one <- function(path, label, kind) {
    f <- read_rds(path)
    pars <- intersect(c("cycle_length", "R", "log10_total0", "b_shape", "sd_iRBC"),
                      f@model_pars)
    out <- map(pars, \(p) grab(f, p)) |> list_rbind() |>
        mutate(fit = label, kind = kind)
    rm(f); invisible(gc())
    out
}

cat("=== reading ===\n")
d <- bind_rows(
    one("_data/wock-schedsim-fit-default-rep1.rds", "sim default-rep1", "simulated"),
    one("_data/wock-schedsim-fit-default-rep2.rds", "sim default-rep2", "simulated"),
    one("_data/wock-schedsim-fit-default-rep3.rds", "sim default-rep3", "simulated"),
    one("_data/wock-fit-np_wide_total0.rds",        "real np_wide_total0", "real"))

cat("\n=== 1. cycle-length bias: does the summary choice move it? ===\n")
cat("Cells: per-trial posterior cycle_length averaged over 13 trials, by\n")
cat("summary statistic, in hours; `bias_*` is that minus the simulated truth\n")
cat("of 45.012 h, positive meaning overestimation. `median_skew` is the\n")
cat("median over trials of each trial's posterior skew: 0 is symmetric.\n")
cat("Simulated fits only -- the real fit has no truth to be biased against.\n\n")
d |>
    filter(par == "cycle_length", kind == "simulated") |>
    summarise(.by = fit,
              mean_cl = mean(mean), median_cl = mean(median),
              bias_mean = mean(mean) - TRUE_CL,
              bias_median = mean(median) - TRUE_CL,
              shift = mean(median) - mean(mean),
              median_skew = stats::median(skew)) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(width = Inf)

cat("\n=== 2. where the summary choice DOES change a number ===\n")
cat("Cells: over each parameter's groups within a fit, the mean of the\n")
cat("per-group posterior means and of the per-group posterior medians, on\n")
cat("the parameter's own scale, and `rel_shift` = (median - mean) / mean.\n")
cat("`median_skew` is the median over groups of the per-group posterior\n")
cat("skew. A parameter with near-zero skew cannot be helped by switching\n")
cat("summary; one with large positive skew has a mean pulled up by a tail.\n\n")
d |>
    summarise(.by = c(fit, par), n = n(),
              mean = mean(mean), median = mean(median),
              rel_shift = (mean(median) - mean(mean)) / mean(mean),
              median_skew = stats::median(skew)) |>
    arrange(par, fit) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(n = Inf, width = Inf)
