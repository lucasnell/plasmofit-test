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

f <- list.files("_data", "^nc-profile-fast-unit[0-9]+[.]rds$", full.names = TRUE)
if (length(f) == 0) stop("no nc-profile-fast-unit*.rds in _data/ -- job unfinished")
r <- bind_rows(lapply(f, readRDS))
cat("units read:", n_distinct(r$unit), "of 14 |", nrow(r), "fits |",
    sprintf("%.1f", sum(r$secs) / 60), "CPU-minutes total\n")
bad <- r |> count(unit, name = "rows") |> filter(rows != max(rows))
if (nrow(bad) > 0) { print(as.data.frame(bad)); stop("units above have the wrong row count") }

verdict <- function(pool, label, knot_name) {
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
cat(sprintf("\nmax |mesh shift| = %.2f log-likelihood units\n", max(abs(C$shift))))
if (max(abs(C$shift)) > 2)
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
