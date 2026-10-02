# Regression test for the package's b_shape change (branch fixed-b-shape).
#
# The change renames the free parameter to `b_shape_free`, makes it zero-sized
# when b_shape is supplied as data, and rebuilds `b_shape` as a transformed
# parameter from whichever source is active. With b_shape = NULL that is
# mathematically a rename: same prior, same bounds, same everything. So a fit
# with b_shape free must describe the same posterior as before the change.
#
# NOT bit-identical draws, for the reason wockner-anchor-regression.R gives at
# length: the references were sampled by a different compiled DSO, and HMC
# amplifies a last-bit difference into a different trajectory within a few
# leapfrog steps. The test that means something is posterior equivalence
# within Monte Carlo error.
#
# Three comparisons:
#   A  np_wide_total0-rebuild vs np_wide_total0   old path, b_shape free
#   B  np_bs400-rebuild       vs np_bs400         old path, tight prior
#   C  np_bs400_data          vs np_bs400         NEW path vs the tight prior
#
# A and B are the regression. C is a different question -- whether fixing
# b_shape as data lands where pinning it with a tight prior does -- and it is
# NOT expected to be exact: sd_log_b_shape = 0.05 is +-10% at 95%, so the
# tight prior is not a point mass.
#
# ON WHAT THIS CAN CONCLUDE. Comparing a shift against its Monte Carlo
# standard error is anti-conservative: two runs also differ in step-size and
# mass-matrix adaptation, which MCSE does not capture. A proper null is two
# runs of IDENTICAL code differing only in sampler seed, and none exists for
# these configs. So large |z| is evidence of a real change, while small |z| is
# consistent with equivalence without establishing it. The script says which
# it found and does not upgrade the second into a PASS.
#
#   srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#       _scripts/wockner-bshape-regression.R 2>&1 | tee _data/bshape-regression.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

PAIRS <- tribble(
    ~cmp, ~new,                      ~old,              ~what,
    "A",  "np_wide_total0-rebuild",  "np_wide_total0",  "old path, b_shape free",
    "B",  "np_bs400-rebuild",        "np_bs400",        "old path, tight prior",
    "C",  "np_bs400_data",           "np_bs400",        "NEW path vs tight prior")

## Posterior mean and its Monte Carlo standard error for every scalar entry.
summ_of <- function(cfg) {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    ss <- rstan::summary(f)$summary
    mp <- f@model_pars
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    health <- list(div = sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
                   max_rhat = max(ss[, "Rhat"], na.rm = TRUE),
                   min_ess = min(ss[, "n_eff"], na.rm = TRUE))
    rm(f); invisible(gc())
    list(tab = tibble(par = rownames(ss), mean = ss[, "mean"],
                      se_mean = ss[, "se_mean"]),
         model_pars = mp, health = health)
}

cat("=== reading fits ===\n")
need <- unique(c(PAIRS$new, PAIRS$old))
S <- map(set_names(need), \(c) { cat("  ", c, "\n"); summ_of(c) })

cat("\n=== 1. structural: parameter name sets ===\n")
cat("Cells: names present in one fit's model_pars and not the other's. The\n")
cat("rename is EXPECTED to add b_shape_free; anything else is a real\n")
cat("structural difference. Not numeric.\n\n")
for (i in seq_len(nrow(PAIRS))) {
    p <- PAIRS[i, ]
    a <- S[[p$new]]$model_pars; b <- S[[p$old]]$model_pars
    cat(sprintf("%s  %s vs %s\n", p$cmp, p$new, p$old))
    cat(sprintf("     only in new: %s\n",
                paste(setdiff(a, b), collapse = ", ") |> (\(x) if (nzchar(x)) x else "(none)")()))
    cat(sprintf("     only in old: %s\n",
                paste(setdiff(b, a), collapse = ", ") |> (\(x) if (nzchar(x)) x else "(none)")()))
}

cat("\n=== 2. posterior shift, scaled by Monte Carlo error ===\n")
cat("Cells: over the scalar entries both fits share, excluding lp__ and any\n")
cat("entry whose se_mean is zero or missing (a parameter held constant has no\n")
cat("Monte Carlo error and no meaningful z). `z` is the difference in\n")
cat("posterior means divided by sqrt(se_mean_new^2 + se_mean_old^2), so it is\n")
cat("the shift in units of the combined Monte Carlo standard error of the two\n")
cat("runs. `n` is how many entries were compared; `max_abs_z` the largest,\n")
cat("`frac_gt2` and `frac_gt3` the share beyond 2 and 3. Dimensionless.\n\n")

z_tab <- map(seq_len(nrow(PAIRS)), \(i) {
    p <- PAIRS[i, ]
    a <- S[[p$new]]$tab; b <- S[[p$old]]$tab
    j <- inner_join(a, b, by = "par", suffix = c("_new", "_old")) |>
        filter(par != "lp__", !startsWith(par, "log_lik"),
               is.finite(se_mean_new), is.finite(se_mean_old),
               se_mean_new > 0, se_mean_old > 0) |>
        mutate(z = (mean_new - mean_old) /
                   sqrt(se_mean_new^2 + se_mean_old^2))
    tibble(cmp = p$cmp, what = p$what, n = nrow(j),
           max_abs_z = max(abs(j$z)), frac_gt2 = mean(abs(j$z) > 2),
           frac_gt3 = mean(abs(j$z) > 3),
           worst = j$par[which.max(abs(j$z))])
}) |> list_rbind()
z_tab |> mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(width = Inf)

cat("\n=== 3. sampler health, both members of each pair ===\n")
cat("Cells: divergent transitions after warmup summed over 4 chains, the\n")
cat("largest split R-hat, and the smallest bulk n_eff, over all parameters.\n")
cat("Absolute, not a comparison.\n\n")
map(need, \(c) tibble(fit = c, div = S[[c]]$health$div,
                      max_rhat = S[[c]]$health$max_rhat,
                      min_ess = S[[c]]$health$min_ess)) |>
    list_rbind() |> mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    print(width = Inf)

cat("\n=== verdict ===\n")
for (i in seq_len(nrow(z_tab))) {
    r <- z_tab[i, ]
    v <- if (r$max_abs_z > 5) {
        "CHANGED -- shift far outside Monte Carlo error"
    } else if (r$max_abs_z > 3) {
        "SUSPECT -- largest shift beyond 3 combined MCSE, look at it"
    } else {
        "consistent with equivalence; NOT a pass (no null run, see header)"
    }
    cat(sprintf("  %s  %-28s %s\n", r$cmp, r$what, v))
}
