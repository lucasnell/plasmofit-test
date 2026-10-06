## Thread 5: read the n_c ladder.
##
## np_wide_both (n_c = 96), np_nc192, np_nc384 differ in n_c and nothing else.
## Raising n_c slows the model's deterministic desynchronisation: the stage
## distribution's sd after the 4.80-cycle Wockner window is 0.224 cycles at
## 96, 0.158 at 192, 0.112 at 384.
##
## The prediction, recorded in wockner-fit.R before these fits were run: if
## b_shape has been absorbing a too-fast decay rate, b_shape FALLS
## monotonically along the ladder. If it is pinned by the data, it stays put.
## Non-monotonic means neither. cycle_length is read alongside, because the
## confound only matters if moving n_c moves the headline number.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-ladder-read.R

suppressPackageStartupMessages({library(rstan); library(dplyr)})

ARMS <- c(np_wide_both = 96L, np_nc192 = 192L, np_nc384 = 384L)

grab <- function(nm) {
    p <- sprintf("_data/wock-fit-%s.rds", nm)
    if (!file.exists(p)) { message("missing: ", p); return(NULL) }
    f <- readRDS(p)
    s <- summary(f)$summary
    pick <- function(par) {
        i <- grep(sprintf("^%s(\\[|$)", par), rownames(s))
        if (length(i) == 0) return(NULL)
        s[i, , drop = FALSE]
    }
    bs <- pick("b_shape"); cl <- pick("cycle_length")
    dv <- sum(rstan::get_divergent_iterations(f))
    tibble(arm = nm, n_c = ARMS[[nm]],
           b_shape_mean = if (is.null(bs)) NA else mean(bs[, "mean"]),
           b_shape_min  = if (is.null(bs)) NA else min(bs[, "mean"]),
           b_shape_max  = if (is.null(bs)) NA else max(bs[, "mean"]),
           cl_mean      = if (is.null(cl)) NA else mean(cl[, "mean"]),
           cl_min       = if (is.null(cl)) NA else min(cl[, "mean"]),
           cl_max       = if (is.null(cl)) NA else max(cl[, "mean"]),
           max_rhat     = max(s[, "Rhat"], na.rm = TRUE),
           min_ess      = min(s[, "n_eff"], na.rm = TRUE),
           div          = dv)
}

res <- bind_rows(lapply(names(ARMS), grab))
cat("=== the n_c ladder ===\n")
cat("Cells: b_shape and cycle_length are posterior means, averaged over the\n",
    "14 (b_shape) and 13 (cycle_length) groups, with the across-group range.\n",
    "cycle_length is in hours. max_rhat/min_ess/div are sampler health over\n",
    "all monitored parameters.\n\n", sep = "")
print(as.data.frame(res), digits = 4)

if (nrow(res) == 3 && !any(is.na(res$b_shape_mean))) {
    b <- res$b_shape_mean; c_ <- res$cl_mean
    cat("\nb_shape     96 -> 192 -> 384:", sprintf("%.2f -> %.2f -> %.2f", b[1], b[2], b[3]), "\n")
    cat("monotone decreasing:", all(diff(b) < 0), "| monotone increasing:", all(diff(b) > 0), "\n")
    cat("cycle_length 96 -> 192 -> 384:", sprintf("%.3f -> %.3f -> %.3f h", c_[1], c_[2], c_[3]), "\n")
    cat("cycle_length total move:", sprintf("%+.3f h", c_[3] - c_[1]), "\n")
}
