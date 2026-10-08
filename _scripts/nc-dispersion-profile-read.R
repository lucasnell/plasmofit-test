## Is the dispersion MAGNITUDE identified, with the initial-spread nuisance
## free? Thread 15; the prerequisite named in claude/ipm-decision.md.
##
## In the Erlang chain, transit coefficient of variation is 1/sqrt(n_c), so
## n_c IS the dispersion magnitude, exactly and monotonically, and b_shape is
## free at every rung. The maximised log-likelihood as a function of n_c is
## therefore a PROFILE LIKELIHOOD in the dispersion with the initial-spread
## nuisance concentrated out. Every rung has the same parameter count, so the
## rungs are directly comparable with no complexity penalty, and units are
## fitted independently so their log-likelihoods add.
##
## PRE-REGISTERED READING RULE. Refined 2026-10-08 BEFORE any of the new
## rungs existed -- a dry run on the three rungs already on disk showed the
## first draft asked the wrong question. It tested for an INTERIOR maximum,
## but a profile that SATURATES toward an asymptote never has one, and
## saturation is the expected shape here. What actually matters is whether
## the 2-log-likelihood interval is BOUNDED in the fine direction, i.e.
## whether the data exclude zero dispersion.
##   TURNS OVER   - argmax is interior and the 2-unit interval excludes both
##                  ends. Dispersion is bounded away from zero AND from large
##                  values: a measured quantity. An IPM would deliver a
##                  parameter with a real interval. Strongest case to build.
##   SATURATES    - the gain at the top rung is < 2 units and the argmax sits
##                  at the fine end. The profile is then consistent with any
##                  dispersion at or below the optimum, INCLUDING ZERO. There
##                  is an upper bound and no lower bound, so an IPM would
##                  return sigma_d against zero with an interval touching it.
##                  That argues AGAINST building one, which is the opposite
##                  of what the structural argument alone suggests.
##   STILL CLIMBING - the gain at the top rung is >= 2 units. The grid is too
##                  coarse to say; extend it before deciding anything.
##   FLAT         - the whole grid spans < 2 units. The magnitude is not
##                  determined either; report the ladder as a sensitivity.
## 2 log-likelihood units because a 95% interval for a single parameter is a
## drop of 1.92; 2 is that rounded, and it is the same currency as everything
## else in thread 15.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-dispersion-profile-read.R | tee _data/nc-dispersion-profile-$(date +%F).txt
##
## The saved output is tracked -- see the exceptions in _data/.gitignore --
## because numbers from it are quoted in the notes.

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

pat <- "^decay-law-(unit|a768-unit|prof-unit)[0-9]+[.]rds$"
f <- list.files("_data", pat, full.names = TRUE)
if (length(f) == 0) stop("no decay-law rds files in _data/")
r <- bind_rows(lapply(f, readRDS)) |> filter(model == "A_sqrt")

## bs_hat was added to the output on 2026-10-08. Files written before that
## have it as NA, which is reported rather than silently dropped.
if (!"bs_hat" %in% names(r)) r$bs_hat <- NA_real_
if (!"cl_hat" %in% names(r)) r$cl_hat <- NA_real_
have_bs <- r |> group_by(n_c) |> summarise(bs = sum(!is.na(bs_hat)), .groups = "drop")

grid <- sort(unique(r$n_c))
cat("units:", n_distinct(r$unit), "| n_c rungs present:",
    paste(grid, collapse = ", "), "\n")
chk <- r |> count(unit, name = "rows") |> filter(rows != length(grid))
if (nrow(chk) > 0) {
    print(as.data.frame(chk))
    stop("unit(s) above are missing rungs -- a job is unfinished or a file is ",
         "stale. Fix before reading, or the pooled profile mixes different ",
         "grids across units.")
}
cat("b_shape recorded at:",
    paste(have_bs$n_c[have_bs$bs > 0], collapse = ", "), "\n\n")

## ---- the pooled profile -----------------------------------------------
pool <- r |> group_by(n_c) |>
    summarise(ll = sum(ll), units = n(), .groups = "drop") |>
    mutate(cv = 1 / sqrt(n_c),
           sd_h = NA_real_,
           d_from_max = ll - max(ll)) |>
    arrange(n_c)

## sd of the stage distribution at the last observation, k = 4.8 cycles, in
## hours, using each unit's own fitted cycle_length averaged over units.
cl_bar <- r |> group_by(n_c) |> summarise(cl = mean(cl_hat, na.rm = TRUE), .groups = "drop")
pool$sd_h <- cl_bar$cl[match(pool$n_c, cl_bar$n_c)] * sqrt(4.8 / pool$n_c)

cat("=== pooled profile likelihood in the dispersion magnitude ===\n")
cat("Cells: ll is the SUM over the 14 grp_init units of model A's maximised\n",
    "profile log-likelihood at that n_c, with b_shape free at every rung and\n",
    "the observation sd profiled out; units are independent so the sum is the\n",
    "joint log-likelihood over all 1130 observations. cv = 1/sqrt(n_c) is the\n",
    "transit coefficient of variation the rung implies. sd_h is the stage\n",
    "distribution's sd at k = 4.8 cycles in hours, at the mean fitted\n",
    "cycle_length. d_from_max is relative to the best rung; 0 marks it.\n\n",
    sep = "")
print(as.data.frame(pool), digits = 5)

best <- pool$n_c[which.max(pool$ll)]
interior <- best != min(grid) && best != max(grid)
within <- pool$n_c[pool$d_from_max > -2]
cat(sprintf("\nbest rung: n_c = %d (cv %.4f) | interior: %s\n",
            best, 1/sqrt(best), interior))
cat(sprintf("within 2 log-likelihood units: n_c in {%s}\n",
            paste(within, collapse = ", ")))
cat(sprintf("profile span over the whole grid: %.1f log-likelihood units\n",
            max(pool$ll) - min(pool$ll)))

bracketed <- interior && min(within) > min(grid) && max(within) < max(grid)

## ---- the b_shape trade-off, which is the actual confound ---------------
cat("\n=== the b_shape / dispersion trade-off ===\n")
bs <- r |> filter(!is.na(bs_hat)) |>
    group_by(n_c) |> summarise(bs_med = median(bs_hat), units = n(),
                               .groups = "drop") |> arrange(n_c)
if (nrow(bs) < 2) {
    cat("Not enough rungs record b_shape yet (added 2026-10-08). Re-run the\n",
        "earlier rungs in profile mode if this is wanted across the grid.\n",
        sep = "")
} else {
    cat("Cells: median over units of the fitted b_shape at each rung. Lower\n",
        "b_shape is a WIDER initial spread. If b_shape falls as n_c rises,\n",
        "the fit is buying back with initial spread what the finer mesh took\n",
        "away in accumulated spread -- that trade-off is the confound.\n\n",
        sep = "")
    print(as.data.frame(bs), digits = 5)
    per <- r |> filter(!is.na(bs_hat)) |> group_by(unit) |>
        summarise(rho = suppressWarnings(cor(log(n_c), log(bs_hat),
                                             method = "spearman")),
                  .groups = "drop")
    cat(sprintf("\nSpearman correlation of log(b_shape) with log(n_c), per unit:\n"))
    cat(sprintf("  median %.2f | negative in %d of %d units\n",
                median(per$rho, na.rm = TRUE), sum(per$rho < 0, na.rm = TRUE),
                sum(!is.na(per$rho))))
}

## ---- verdict ------------------------------------------------------------
cat("\n=== pre-registered verdict ===\n")
top <- pool$n_c[nrow(pool)]
top_gain <- pool$ll[nrow(pool)] - pool$ll[nrow(pool) - 1]
cat(sprintf("gain at the top rung (%d over %d): %+.2f log-likelihood units\n",
            top, pool$n_c[nrow(pool) - 1], top_gain))
if (max(pool$ll) - min(pool$ll) < 2) {
    cat("FLAT. The profile spans less than 2 log-likelihood units over the\n",
        "whole grid, so the dispersion magnitude is not determined either.\n",
        "Neither a free stage-rate parameter nor an IPM helps. Report the\n",
        "ladder as a sensitivity analysis and say the rate is unidentified.\n",
        sep = "")
} else if (bracketed) {
    cat("TURNS OVER. The pooled profile has an interior maximum and the\n",
        "2-log-likelihood interval excludes both ends of the grid, so the\n",
        "dispersion rate is bounded away from zero and from large values: a\n",
        "MEASURED quantity. An IPM would deliver a parameter with a real\n",
        "interval, and a free stage-rate parameter reaches it too as long as\n",
        "the optimum sits ABOVE the Erlang floor at the chosen mesh.\n",
        sep = "")
} else if (top_gain < 2) {
    cat("SATURATES. The argmax is at the fine end of the grid and the top\n",
        "rung buys less than 2 log-likelihood units, so the profile is\n",
        "approaching an asymptote. The data put an UPPER bound on dispersion\n",
        "and no lower bound -- they do not exclude zero. An IPM would return\n",
        "sigma_d pressed against zero with an interval touching it, which is\n",
        "a boundary estimate rather than a rate. This argues AGAINST building\n",
        "one, opposite to what the structural argument alone suggests, and\n",
        "for picking n_c by elpd and reporting the ladder.\n", sep = "")
} else {
    cat("STILL CLIMBING. The top rung still buys 2 or more log-likelihood\n",
        "units, so the grid is too coarse to say where the optimum is or\n",
        "whether there is one. Extend the grid before deciding anything.\n",
        sep = "")
}
