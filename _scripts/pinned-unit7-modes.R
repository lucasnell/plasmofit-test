## Follow-up to _scripts/pinned-validate.R: the pinned refit of entry 43 failed
## the gate only through unit 7 (MMV048_PIB). Does chain 1 sit in a second
## mode, or move between modes?
##
## Cells: per chain, over its 1000 post-warmup draws of cycle_length[7] (h):
## mean, the share of draws below 42.6 h (between the two clusters seen in
## the per-chain means), and that share in each quarter of the chain, in
## draw order. lp__ by chain for the same fit. Fits: SLURM 30830_43
## (original) and 30856 (pinned).
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/pinned-unit7-modes.R > _data/pinned-unit7-modes-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))
CUT <- 42.6
for (nm in c("np_bs15_nc384", "np_bs15_nc384-pinned")) {
    f <- readRDS(sprintf("_data/wock-fit-%s.rds", nm))
    a <- extract(f, "cycle_length[7]", permuted = FALSE)[, , 1]
    lp <- extract(f, "lp__", permuted = FALSE)[, , 1]
    cat("==", nm, "\n")
    for (ch in seq_len(ncol(a))) {
        x <- a[, ch]; q <- split(x, rep(1:4, each = length(x) / 4))
        cat(sprintf("chain %d: mean %.3f | below %.1f h: %5.1f%% | by quarter %s | lp__ %.1f\n",
                    ch, mean(x), CUT, 100 * mean(x < CUT),
                    paste(sprintf("%3.0f%%", sapply(q, function(v) 100 * mean(v < CUT))), collapse = " "),
                    mean(lp[, ch])))
    }
    cat("  per-chain quantiles 5/50/95% and mean by quarter (h):\n")
    for (ch in seq_len(ncol(a))) {
        x <- a[, ch]; q <- split(x, rep(1:4, each = length(x) / 4))
        cat(sprintf("  chain %d: %s | quarters %s\n", ch,
                    paste(sprintf("%.2f", quantile(x, c(.05, .5, .95))), collapse = "/"),
                    paste(sprintf("%.2f", sapply(q, mean)), collapse = " ")))
    }
    lo <- a < CUT
    if (any(lo)) cat(sprintf("mean lp__ when below %.1f: %.1f; above: %.1f\n", CUT,
                             mean(lp[lo]), mean(lp[!lo])))
    cat("\n")
}
