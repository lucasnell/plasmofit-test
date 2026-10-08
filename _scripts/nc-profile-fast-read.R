## Read the dense dispersion profile. Thread 15; companion to
## _scripts/nc-profile-fast.R and the decision in claude/ipm-decision.md.
##
## THE READING RULE IS THE ONE PRE-REGISTERED IN
## _scripts/nc-dispersion-profile-read.R, before any of this existed:
##   TURNS OVER     - interior maximum and the 2-log-likelihood interval
##                    excludes both ends. Dispersion is bounded away from zero
##                    and from large values: a MEASURED quantity.
##   SATURATES      - argmax at the fine end and the top rung buys < 2 units.
##                    Upper bound only; zero is not excluded, so an IPM would
##                    return a boundary estimate rather than a rate.
##   STILL CLIMBING - the top rung still buys >= 2 units; the grid is too
##                    coarse to say.
##   FLAT           - the whole grid spans < 2 units.
## 2 units because a 95% interval for one parameter is a drop of 1.92.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-profile-fast-read.R | tee _data/nc-profile-fast-$(date +%F).txt
##
## The saved output is tracked -- see the exceptions in _data/.gitignore.

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

## Two runs may be on disk: the original two-start run (untagged) and the
## eight-start re-run (-s8). The eight-start set is the one to read, because
## more starts can only raise a maximised log-likelihood. The DIFFERENCE
## between them is a direct measurement of how far the two-start optimiser
## stopped short, which is worth reporting in its own right.
read_set <- function(pat) {
    fs <- list.files("_data", pat, full.names = TRUE)
    if (length(fs) == 0) return(NULL)
    bind_rows(lapply(fs, readRDS))
}
r2 <- read_set("^nc-profile-fast-unit[0-9]+[.]rds$")
r8 <- read_set("^nc-profile-fast-unit[0-9]+-s8[.]rds$")
if (is.null(r2) && is.null(r8)) stop("no nc-profile-fast-unit*.rds in _data/")
r <- if (!is.null(r8)) r8 else r2
which_set <- if (!is.null(r8)) "eight starts" else "two starts"
cat("reading the", which_set, "run |", n_distinct(r$unit), "of 14 units |",
    nrow(r), "fits |", sprintf("%.1f", sum(r$secs)/60), "CPU-minutes\n")
if (!is.null(r8) && !is.null(r2) && n_distinct(r8$unit) == n_distinct(r2$unit)) {
    cmp <- inner_join(r2 |> select(unit, model, knot, ll2 = ll),
                      r8 |> select(unit, model, knot, ll8 = ll),
                      by = c("unit", "model", "knot")) |>
           mutate(gain = ll8 - ll2)
    g <- cmp |> group_by(model, knot) |>
         summarise(pooled_gain = sum(gain), worst_unit = max(gain), .groups = "drop") |>
         arrange(desc(pooled_gain))
    cat("\n=== what the extra six starts bought ===\n")
    cat("Cells: pooled_gain is the sum over units of (eight-start ll minus\n",
        "two-start ll) at that rung, which can only be >= 0; worst_unit is the\n",
        "largest single-unit gain there. Anything near the 2-unit currency means\n",
        "the two-start profile was shaped by where the optimiser stopped.\n\n", sep = "")
    print(as.data.frame(head(g, 10)), digits = 4)
    cat(sprintf("\ntotal gain over all rungs: %.1f ll units | rungs gaining > 2: %d of %d\n",
                sum(cmp$gain), sum(g$pooled_gain > 2), nrow(g)))
    if (any(cmp$gain < -1e-6))
        cat("WARNING:", sum(cmp$gain < -1e-6), "cells went DOWN with more starts,",
            "which is impossible and means the two sets are not comparable.\n")
}
bad <- r |> count(unit, name = "rows") |> filter(rows != max(rows))
if (nrow(bad) > 0) { print(as.data.frame(bad)); stop("units above have the wrong row count") }

## A profile likelihood in a smooth parameter should be SMOOTH, so deviation
## from a smooth fit measures how far the optimiser stopped short. This has to
## be checked BEFORE any verdict, because the rule works in 2-log-likelihood
## units and under-optimisation of that size would invent or erase a turnover.
## Learned the hard way twice now: the decay-law test (29684) had the same
## problem, where six nested models scored worse than models they contain.
roughness <- function(pool, knot_name) {
    k <- pool[[knot_name]]
    if (length(k) < 6) return(NA_real_)
    r <- pool$ll - predict(loess(ll ~ log(k), data = data.frame(ll = pool$ll, k = k),
                                 span = 0.75, degree = 2))
    max(abs(r))
}

verdict <- function(pool, label, knot_name) {
    rough <- roughness(pool, knot_name)
    best <- pool[[knot_name]][which.max(pool$ll)]
    grid <- pool[[knot_name]]
    within <- grid[pool$ll - max(pool$ll) > -2]
    interior <- best != min(grid) && best != max(grid)
    bracketed <- interior && min(within) > min(grid) && max(within) < max(grid)
    top_gain <- pool$ll[nrow(pool)] - pool$ll[nrow(pool) - 1]
    cat(sprintf("\n-- %s --\n", label))
    cat(sprintf("best %s = %g | 2-unit interval {%g, %g} | span %.1f | top-rung gain %+.2f\n",
                knot_name, best, min(within), max(within),
                max(pool$ll) - min(pool$ll), top_gain))
    cat(sprintf("roughness (max |deviation from a smooth profile|) = %.2f ll units\n", rough))
    if (!is.na(rough) && rough > 2) {
        cat("NO VERDICT. The profile is rougher than the 2-log-likelihood\n",
            "currency the rule is written in, so the interval above is not\n",
            "trustworthy and a turnover could be invented or erased by where\n",
            "the optimiser stopped. Under-optimisation can only push ll DOWN,\n",
            "so rungs with large negative residuals are the suspects. Refit\n",
            "with more starts (see _scripts/profile-noise-check.R) and re-read.\n",
            sep = "")
        return(invisible(NULL))
    }
    if (max(pool$ll) - min(pool$ll) < 2) {
        cat("FLAT: the magnitude is not determined. Neither fix helps.\n")
    } else if (bracketed) {
        cat("TURNS OVER: dispersion is bounded away from zero AND from large\n",
            "values, so it is a measured quantity. An IPM would deliver a\n",
            "parameter with a real interval -- the strongest case to build one.\n",
            sep = "")
    } else if (top_gain < 2) {
        cat("SATURATES: upper bound on dispersion, no lower bound, so the data\n",
            "do not exclude zero. An IPM would return a boundary estimate, not\n",
            "a rate. Argues AGAINST building one.\n", sep = "")
    } else {
        cat("STILL CLIMBING: the grid is too coarse to say. Extend it.\n")
    }
}

cat("\n=== A. the exact chain, nine integer rungs ===\n")
cat("Cells: ll is the SUM over the 14 grp_init units of the maximised profile\n",
    "log-likelihood with b_shape free and the observation sd profiled out;\n",
    "units are independent so the sum is the joint log-likelihood over all\n",
    "1130 observations. cv = 1/sqrt(n_c) is the implied transit coefficient of\n",
    "variation. Every rung has the same parameter count. No priors.\n\n", sep = "")
A <- r |> filter(model == "chain") |> group_by(n_c = knot) |>
    summarise(ll = sum(ll), cl_bar = mean(cl_hat), bs_med = median(bs_hat),
              .groups = "drop") |>
    mutate(cv = 1 / sqrt(n_c), d_max = ll - max(ll)) |> arrange(n_c)
print(as.data.frame(A), digits = 5)
verdict(A, "chain", "n_c")

cat("\n\n=== B. the gamma IPM: mesh FIXED at 192, n_eff continuous ===\n")
cat("Cells: as above, but the mesh is held fixed and n_eff -- the chain's\n",
    "dispersion parameter, made continuous -- is varied. This is the forward\n",
    "map an IPM would use. If the chain profile and this one disagree about\n",
    "where the optimum is, the chain's answer was partly a mesh effect.\n\n",
    sep = "")
B <- r |> filter(model == "ipm_gamma") |> group_by(n_eff = knot) |>
    summarise(ll = sum(ll), cl_bar = mean(cl_hat), bs_med = median(bs_hat),
              .groups = "drop") |>
    mutate(cv = 1 / sqrt(n_eff), d_max = ll - max(ll)) |> arrange(n_eff)
print(as.data.frame(B), digits = 5)
verdict(B, "gamma IPM at fixed mesh", "n_eff")

cat("\n\n=== C. is the fixed mesh fine enough? ===\n")
cat("Cells: the same n_eff fitted at mesh 192 and at mesh 384, summed over\n",
    "units. A large shift would mean the rows above are mesh artefacts rather\n",
    "than results. Judge against the 2-log-likelihood currency used throughout.\n\n",
    sep = "")
C <- r |> filter(model %in% c("ipm_gamma", "ipm_gamma_meshcheck")) |>
    group_by(n_eff = knot, M) |> summarise(ll = sum(ll), .groups = "drop") |>
    pivot_wider(names_from = M, values_from = ll, names_prefix = "M") |>
    filter(!is.na(M384)) |> mutate(shift = M384 - M192)
print(as.data.frame(C), digits = 6)
valid <- C[!is.na(C$shift), , drop = FALSE]
if (nrow(valid) == 0) {
    cat("\nNo rung was run at BOTH meshes, so the mesh check is empty. The\n",
        "n_eff values in NE_CHECK must be drawn from NE_IPM for the two to\n",
        "pair up; they were not. Fix and re-run before trusting section B.\n",
        sep = "")
} else {
    cat(sprintf("\nmax |mesh shift| = %.2f ll units over %d rung(s) run at both meshes\n",
                max(abs(valid$shift)), nrow(valid)))
}
if (nrow(valid) > 0 && max(abs(valid$shift)) > 2)
    cat("WARNING: larger than the 2-unit currency. Treat section B as mesh-dependent.\n")

cat("\n\n=== D. the cycle_length that comes with each rung ===\n")
cat("Cells: mean over the 14 units of the fitted cycle_length in hours, and\n",
    "the median fitted b_shape, at each rung of both profiles. This is the\n",
    "quantity the project is actually trying to report, so its spread across\n",
    "the 2-log-likelihood interval is the honest uncertainty from this\n",
    "assumption alone.\n\n", sep = "")
D <- bind_rows(A |> transmute(model = "chain", knot = n_c, ll, cl_bar, bs_med),
               B |> transmute(model = "ipm_gamma", knot = n_eff, ll, cl_bar, bs_med))
for (m in unique(D$model)) {
    s <- D |> filter(model == m); w <- s |> filter(ll - max(s$ll) > -2)
    cat(sprintf("%-10s : cycle_length over the 2-unit interval = %.2f to %.2f h (spread %.2f)\n",
                m, min(w$cl_bar), max(w$cl_bar), diff(range(w$cl_bar))))
}
print(as.data.frame(D), digits = 5)
