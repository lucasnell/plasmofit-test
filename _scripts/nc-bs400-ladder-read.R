## Read the production n_c ladder with b_shape pinned at 400: entries 39 and
## 40, both reseeded (SLURM 30525 and 30829, WOCKFIT_SUFFIX=-seed2).
##
## np_bs400_nc192 and np_bs400_nc384 differ in n_c and nothing else. The
## question: does the n_c effect on cycle_length survive a pinned b_shape?
##
## The prediction, recorded in claude/findings.md before 30829 ran ("A
## prediction for the retry"): the effect PERSISTS and is larger at b_shape
## 400 than at 15. The least-squares forward-map table at b_shape 400 gives a
## 192 -> 384 shift of -1.21 h (48 h window), -0.94 (96 h), -0.72 (144 h),
## -0.95 (96 h, nuisances free). If the run converges and the shift is near
## zero, the mechanistic account is wrong.
##
## Order of reading, per the notes: lp__ per chain first, then the R-hat
## gate (max < 1.05), then cycle_length. b_shape is data here, so its R-hat is
## NA and is excluded; the count of NA R-hats is printed so that is visible.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-bs400-ladder-read.R \
##       > _data/nc-bs400-ladder-read-2026-10-10.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))

RUNGS <- c(nc192 = "np_bs400_nc192-seed2", nc384 = "np_bs400_nc384-seed2")
GATE <- 1.05

fits <- lapply(RUNGS, function(nm) readRDS(sprintf("_data/wock-fit-%s.rds", nm)))
dats <- lapply(RUNGS, function(nm) readRDS(sprintf("_data/wock-data-%s.rds", nm)))

## Entry-wise comparison is only valid if both fits derived the same group
## levels in the same order (claude/gotchas.md).
stopifnot(identical(attr(dats$nc192, "levels"), attr(dats$nc384, "levels")),
          dats$nc192$n_c == 192L, dats$nc384$n_c == 384L,
          all(dats$nc192$b_shape_data == 400), all(dats$nc384$b_shape_data == 400),
          identical(dats$nc192$y, dats$nc384$y))
cl_levels <- attr(dats$nc192, "levels")[["grp_cl"]]
if (is.null(cl_levels)) cl_levels <- attr(dats$nc192, "levels")[[3]]

health <- function(f, nm) {
    s <- summary(f)$summary
    lp <- extract(f, "lp__", permuted = FALSE)[, , 1, drop = FALSE]
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    mtd <- f@stan_args[[1]]$control$max_treedepth
    if (is.null(mtd)) mtd <- 10L
    cat(sprintf("--- %s ---\n", nm))
    cat(sprintf("chains %d x %d post-warmup draws\n", dim(lp)[2], dim(lp)[1]))
    cat("lp__ mean by chain: ",
        paste(sprintf("%.1f", apply(lp, 2, mean)), collapse = ", "),
        sprintf("  (spread %.1f)\n", diff(range(apply(lp, 2, mean)))))
    cat("divergent by chain: ",
        paste(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)), collapse = ", "),
        sprintf("  (total %d of %d, %.1f%%)\n",
                sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1))),
                sum(vapply(sp, nrow, integer(1))),
                100 * mean(unlist(lapply(sp, function(x) x[, "divergent__"])))))
    td <- unlist(lapply(sp, function(x) x[, "treedepth__"]))
    cat(sprintf("treedepth at max (%d): %.1f%% of transitions\n", mtd, 100 * mean(td >= mtd)))
    rh <- s[, "Rhat"]
    cat(sprintf("max R-hat %.4f over %d quantities (%d NA, excluded) -> %s\n",
                max(rh, na.rm = TRUE), sum(!is.na(rh)), sum(is.na(rh)),
                if (max(rh, na.rm = TRUE) < GATE) "PASS" else "FAIL"))
    worst <- head(sort(rh, decreasing = TRUE), 5)
    cat("five worst R-hat: ", paste(sprintf("%s %.3f", names(worst), worst), collapse = "; "), "\n")
    cat(sprintf("min n_eff %.0f\n\n", min(s[, "n_eff"], na.rm = TRUE)))
    invisible(max(rh, na.rm = TRUE))
}

cat("=== sampler health, read lp__ per chain before the gate ===\n\n")
rhat <- mapply(health, fits, RUNGS)

## Per-trial cycle_length.
cl_draws <- lapply(fits, function(f) extract(f, "cycle_length")[[1]])   # draws x 13
cl_mean  <- sapply(cl_draws, colMeans)
cl_sd    <- sapply(cl_draws, function(x) apply(x, 2, sd))
cl_rhat  <- sapply(fits, function(f) {
    s <- summary(f, pars = "cycle_length")$summary; s[, "Rhat"] })
stopifnot(nrow(cl_mean) == 13L)

cat("=== cycle_length by trial ===\n")
cat("Cells: posterior mean (and posterior sd) of cycle_length in hours for each\n",
    "of the 13 trials (grp_cl), b_shape pinned at 400. diff = nc384 - nc192 for\n",
    "the SAME trial (paired by level name; levels verified identical); negative\n",
    "means the finer rung gives a shorter cycle, the predicted direction.\n\n", sep = "")
tab <- data.frame(trial = if (length(cl_levels) == 13) cl_levels else seq_len(13),
                  nc192 = sprintf("%.3f (%.3f)", cl_mean[, 1], cl_sd[, 1]),
                  nc384 = sprintf("%.3f (%.3f)", cl_mean[, 2], cl_sd[, 2]),
                  diff  = sprintf("%+.3f", cl_mean[, 2] - cl_mean[, 1]),
                  rhat192 = sprintf("%.3f", cl_rhat[, 1]),
                  rhat384 = sprintf("%.3f", cl_rhat[, 2]))
print(tab, row.names = FALSE)

## Mean over trials, computed per draw so a posterior sd exists for it.
tm <- lapply(cl_draws, rowMeans)
d  <- cl_mean[, 2] - cl_mean[, 1]
cat("\n=== summary over the 13 trials ===\n")
cat("Cells: hours. 'mean over trials' is the posterior mean of the per-draw\n",
    "mean over 13 trials, with its posterior sd. The two fits are separate runs\n",
    "on the SAME data, so the difference is a model-assumption effect, not a\n",
    "sampling comparison, and no sd is attached to it.\n\n", sep = "")
cat(sprintf("mean over trials, n_c 192: %.3f (posterior sd %.3f), range %.3f to %.3f\n",
            mean(tm$nc192), sd(tm$nc192), min(cl_mean[, 1]), max(cl_mean[, 1])))
cat(sprintf("mean over trials, n_c 384: %.3f (posterior sd %.3f), range %.3f to %.3f\n",
            mean(tm$nc384), sd(tm$nc384), min(cl_mean[, 2]), max(cl_mean[, 2])))
cat(sprintf("paired shift 192 -> 384: mean %+.3f h, median %+.3f, range %+.3f to %+.3f\n",
            mean(d), median(d), min(d), max(d)))
cat(sprintf("trials shifting down: %d of 13\n", sum(d < 0)))
cat("forward-map prediction at b_shape 400 (findings.md): -0.72 to -1.21 h, -0.94 at the 96 h window\n")

## Per-chain cycle_length, so a chain that converged to a different mode is
## visible even when R-hat passes.
cat("\n=== mean over trials of cycle_length, by chain ===\n")
cat("Cells: per-chain posterior mean of the mean over 13 trials, hours.\n")
for (k in names(fits)) {
    a <- extract(fits[[k]], "cycle_length", permuted = FALSE)
    cat(sprintf("%s: %s\n", k, paste(sprintf("%.3f", apply(a, 2, mean)), collapse = ", ")))
}

## For context only, the non-converged originals from 29635, which are not
## posteriors and are not reported: is the lp__ level of the reseed where the
## three agreeing chains were, or where the stuck one was?
cat("\n=== context: 29635 originals, NOT posteriors ===\n")
for (nm in c("np_bs400_nc192", "np_bs400_nc384")) {
    p <- sprintf("_data/wock-fit-%s.rds", nm)
    if (!file.exists(p)) { cat("missing:", p, "\n"); next }
    lp <- extract(readRDS(p), "lp__", permuted = FALSE)[, , 1, drop = FALSE]
    cat(sprintf("%s lp__ by chain: %s\n", nm,
                paste(sprintf("%.1f", apply(lp, 2, mean)), collapse = ", ")))
}
