## Planning aid (2026-10-27) for a single-unit study with b_shape always
## fixed, at the truth or above or below it, crossed with fit n_c at, above,
## or below the generating n_c. Not a study script.
##
## 1. Spread matching. Total stage-age variance (cycles^2) at time t under
##    the model: var_init(b) + k / n_c, with var_init(b) = 1 / (4 (2b + 1))
##    and k = t / cycle_length (unitsim.R header). For the main truth (gen
##    n_c 384, b_shape 20.10, cycle_length 42.844 h), the b_shape a fit at
##    n_c' needs to reproduce the true total spread at t = 168 and 192 h
##    (NA = impossible: the fit's accumulated spread alone exceeds the
##    truth's total), and the initial spread sd (cycles) of candidate pins.
## 2. Cost: hours per chain (warmup + sampling) of the main study's fixed
##    fits at fit n_c 384, 2000 + 2000 iterations.
##
##   Rscript --vanilla _scripts/unitsim-pinplan.R > _data/unitsim-read-pinplan-<date>.txt
suppressPackageStartupMessages(library(dplyr))
vi <- function(b) 1 / (4 * (2 * b + 1)); bv <- function(v) ifelse(v > 0, (1 / (4 * v) - 1) / 2, NA)
b0 <- 20.098; cl <- 42.844; g <- 384
cat("1. Spread-matched b_shape by fit n_c (truth: gen n_c 384, b_shape 20.10)\n\n")
cat(" fit_nc  sd_total_168  b_match_168  b_match_192\n")
for (nc in c(192, 256, 288, 320, 384, 512, 768, 1024)) {
    tot <- function(t) vi(b0) + (t / cl) / g
    cat(sprintf(" %6d  %12.4f  %11.2f  %11.2f\n", nc, sqrt(tot(168)),
                bv(tot(168) - (168 / cl) / nc), bv(tot(192) - (192 / cl) / nc)))
}
cat("\nInitial spread sd (cycles) of candidate pins; sd of the truth's total spread at 168 h = ",
    sprintf("%.4f", sqrt(vi(b0) + (168 / cl) / g)), "\n\n", sep = "")
for (b in c(5.86, 10.05, 15, 20.098, 40.2, 80.4, 400))
    cat(sprintf(" b_shape %7.2f  initial sd %.4f  total sd at 168 h, fit n_c 384: %.4f\n",
                b, sqrt(vi(b)), sqrt(vi(b) + (168 / cl) / 384)))
cat("\n2. Hours per chain, main study's fixed fits (fit n_c 384, 2000 + 2000 iterations)\n\n")
s <- bind_rows(lapply(Sys.glob("_data/unitsim/study/sum-w*-fix-rep*.rds"), readRDS)) |>
    filter(is.na(prior) | prior == "prod", is.na(truth_set) | truth_set == "main")
h <- unlist(lapply(s$per_chain, \(p) p$hours))
cat(sprintf(" %d fits, %d chains: hours per chain median %.3f, 90%% %.3f, max %.3f; ms per gradient median %.2f\n",
            nrow(s), length(h), median(h), quantile(h, 0.9), max(h),
            median(unlist(lapply(s$per_chain, \(p) p$ms_per_grad)))))
