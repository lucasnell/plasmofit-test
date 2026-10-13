## Read the b_shape ladder at n_c = 384: SLURM 30830 (entries 42-48), 30837
## (49-53) and 30842 (54), plus entry 40's reseed (30829) as the failed
## default-settings 400 rung. Reads whatever has landed; rerun as more do.
##
## The rung choice and prediction were recorded in _scripts/wockner-fit.R
## above entries 43-48 before any fit ran:
##   - lower rungs CONVERGE under the settings that failed at 400, and the
##     elpd plateau sits BELOW 400 -> the pin is what made 384 hard, and the
##     400 pin is conditional on n_c = 96;
##   - rungs fail regardless of b_shape -> the difficulty is n_c-driven.
## cycle_length is reported as its spread across CONVERGED rungs, never the
## value at the best one.
##
## Order of reading: lp__ by chain, then the gate (max R-hat < 1.05; b_shape
## is data in pinned rungs, so its NA R-hats are excluded and counted). Only
## fits that pass enter the elpd and cycle_length sections.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc384-ladder-read.R \
##       > _data/nc384-ladder-read-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(loo)})

GATE <- 1.05
## name -> (b_shape rung, settings); "free" is estimated b_shape.
FITS <- list(
    `np_bs400_nc384-seed2`    = c("400",  "default"),
    np_bs15_nc384             = c("15",   "default"),
    np_bs50_nc384             = c("50",   "default"),
    np_bs100_nc384            = c("100",  "default"),
    np_bs250_nc384            = c("250",  "default"),
    np_bs600_nc384            = c("600",  "default"),
    np_wide_total0_nc384      = c("free", "default"),
    np_bs400_nc384_t          = c("400",  "tuned"),
    np_bs15_nc384_t           = c("15",   "tuned"),
    np_bs50_nc384_t           = c("50",   "tuned"),
    np_bs100_nc384_t          = c("100",  "tuned"),
    np_bs250_nc384_t          = c("250",  "tuned"),
    np_bs600_nc384_t          = c("600",  "tuned"),
    np_wide_total0_nc384_t    = c("free", "tuned"))

path <- function(kind, nm) sprintf("_data/wock-%s-%s.rds", kind, nm)
have <- names(FITS)[vapply(names(FITS), function(nm)
    all(file.exists(path(c("fit", "data", "fit-LOO"), nm))), logical(1))]
cat("fits present:", length(have), "of", length(FITS), "\n")
cat("missing:", paste(setdiff(names(FITS), have), collapse = ", "), "\n\n")

## Entry-wise comparison needs identical group levels and data.
ref <- readRDS(path("data", have[1]))
for (nm in have) {
    d <- readRDS(path("data", nm))
    stopifnot(identical(attr(d, "levels"), attr(ref, "levels")),
              identical(d$y, ref$y), d$n_c == 384L)
    if (FITS[[nm]][1] != "free")
        stopifnot(d$b_shape_fixed == 1L, all(d$b_shape_data == as.numeric(FITS[[nm]][1])))
    else stopifnot(d$b_shape_fixed == 0L)
}
cl_levels <- attr(ref, "levels")[["trial"]]   # grp_cl = trial

cat("=== sampler health: read lp__ by chain before the gate ===\n")
cat("Cells: per-chain posterior mean of lp__ (4 chains x 1000 post-warmup\n",
    "draws); spread = max - min of those; max R-hat over monitored quantities\n",
    "excluding NA (b_shape when pinned); divergent and treedepth-at-cap as % of\n",
    "4000 post-warmup transitions; hours = longest chain, warmup + sampling.\n\n", sep = "")
H <- list(); fits <- list()
for (nm in have) {
    f <- readRDS(path("fit", nm)); fits[[nm]] <- f
    s <- summary(f)$summary
    lp <- extract(f, "lp__", permuted = FALSE)[, , 1, drop = FALSE]
    lpm <- apply(lp, 2, mean)
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    mtd <- f@stan_args[[1]]$control$max_treedepth; if (is.null(mtd)) mtd <- 10L
    rh <- s[, "Rhat"]; mx <- max(rh, na.rm = TRUE)
    H[[nm]] <- data.frame(
        fit = nm, b_shape = FITS[[nm]][1], settings = FITS[[nm]][2],
        lp_by_chain = paste(sprintf("%.1f", lpm), collapse = ", "),
        lp_spread = round(diff(range(lpm)), 1),
        max_rhat = round(mx, 4), n_na = sum(is.na(rh)),
        gate = if (mx < GATE) "PASS" else "FAIL",
        div_pct = round(100 * mean(unlist(lapply(sp, function(x) x[, "divergent__"]))), 2),
        td_cap_pct = round(100 * mean(unlist(lapply(sp, function(x) x[, "treedepth__"])) >= mtd), 1),
        min_neff = round(min(s[, "n_eff"], na.rm = TRUE)),
        hours = round(max(rowSums(get_elapsed_time(f))) / 3600, 1))
}
H <- do.call(rbind, H); rownames(H) <- NULL
print(H[, c("fit", "lp_by_chain", "lp_spread", "max_rhat", "gate")], row.names = FALSE)
cat("\n")
print(H[, c("fit", "b_shape", "settings", "n_na", "div_pct", "td_cap_pct", "min_neff", "hours")],
      row.names = FALSE)

worst <- function(f) { rh <- summary(f)$summary[, "Rhat"]
    w <- head(sort(rh, decreasing = TRUE), 3)
    paste(sprintf("%s %.3f", names(w), w), collapse = "; ") }
cat("\nworst R-hats, failing fits only:\n")
for (nm in H$fit[H$gate == "FAIL"]) cat(sprintf("  %s: %s\n", nm, worst(fits[[nm]])))

ok <- H$fit[H$gate == "PASS"]
cat(sprintf("\nconverged: %d of %d present: %s\n", length(ok), nrow(H),
            paste(ok, collapse = ", ")))

## Free-b_shape rungs: where does b_shape go at n_c = 384?
free <- intersect(H$fit[H$b_shape == "free"], names(fits))
if (length(free)) {
    cat("\n=== b_shape where it is estimated ===\n")
    cat("Cells: posterior mean of b_shape per grp_init unit (14), summarised\n",
        "across units. Only meaningful for a fit that passes the gate.\n", sep = "")
    for (nm in free) {
        b <- colMeans(extract(fits[[nm]], "b_shape")[[1]])
        cat(sprintf("  %s [%s]: mean %.2f, median %.2f, range %.2f to %.2f\n",
                    nm, H$gate[H$fit == nm], mean(b), median(b), min(b), max(b)))
    }
}

if (length(ok) == 0) { cat("\nNo converged fit; nothing further is read.\n"); quit(save = "no") }

## cycle_length across converged rungs.
cat("\n=== cycle_length, converged fits only ===\n")
cat("Cells: hours. 'mean' is the posterior mean of the per-draw mean over 13\n",
    "trials (grp_cl), 'sd' its posterior sd; range is across the 13 per-trial\n",
    "posterior means. Same data in every fit, so differences between rows are\n",
    "model-assumption effects, not sampling comparisons.\n\n", sep = "")
CL <- do.call(rbind, lapply(ok, function(nm) {
    x <- extract(fits[[nm]], "cycle_length")[[1]]; tm <- rowMeans(x); pm <- colMeans(x)
    data.frame(fit = nm, b_shape = FITS[[nm]][1], mean = round(mean(tm), 3),
               sd = round(sd(tm), 3), min_trial = round(min(pm), 3),
               max_trial = round(max(pm), 3))
}))
print(CL, row.names = FALSE)
cat(sprintf("\nspread of the mean over converged rungs: %.3f to %.3f h (%.3f h)\n",
            min(CL$mean), max(CL$mean), diff(range(CL$mean))))

per_trial <- sapply(ok, function(nm) colMeans(extract(fits[[nm]], "cycle_length")[[1]]))
rownames(per_trial) <- cl_levels
cat("\nper-trial posterior means of cycle_length (h), converged fits:\n")
print(round(per_trial, 2))

## Paired PSIS-LOO over the same 1130 observations.
if (length(ok) >= 2) {
    L <- lapply(ok, function(nm) readRDS(path("fit-LOO", nm))$observation)
    names(L) <- ok
    cat("\n=== paired PSIS-LOO, observation level, converged fits only ===\n")
    cat("Cells: elpd_diff is elpd of the row minus elpd of the best fit, summed\n",
        "over the same 1130 observations (so 0 = best, more negative = worse),\n",
        "with the standard error of the paired difference. Pareto k > 0.7\n",
        "counts are per fit.\n\n", sep = "")
    print(loo_compare(L), simplify = FALSE, digits = 2)
    cat("\nPareto k > 0.7:\n")
    for (nm in ok) cat(sprintf("  %s: %d\n", nm, sum(pareto_k_values(L[[nm]]) > 0.7)))
    fr <- intersect(ok, H$fit[H$b_shape == "free"])
    if (length(fr)) {
        cat("\npaired against the free rung (", fr[1], "); positive = pinned rung better:\n", sep = "")
        for (nm in setdiff(ok, fr)) {
            d <- L[[nm]]$pointwise[, "elpd_loo"] - L[[fr[1]]]$pointwise[, "elpd_loo"]
            cat(sprintf("  %-24s %+8.2f (se %.2f)\n", nm, sum(d), sqrt(length(d)) * sd(d)))
        }
    }
}
