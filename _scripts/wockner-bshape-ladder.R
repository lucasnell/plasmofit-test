# The b_shape ladder: what does fixing synchrony cost, and what does it move?
#
# b_shape is the shape of a symmetric Beta over cycle position, so it is a
# synchrony parameter and converts to a starting-stage age range (see
# findings.md, "What b_shape means biologically"). It is unidentified over at
# least 15-65 and its effect collapses with its value, so the question is not
# what it is but what fixing it costs.
#
# Five rungs, b_shape pinned by a tight prior (sd_log_b_shape = 0.05) at 50,
# 84, 100, 150, 250, all carrying the corrected log10_total0 prior. They are
# read against np_wide_total0, which has the same log10_total0 prior and a
# FREE b_shape, and against np_wide_both, which also widens b_shape.
#
# max_shape is 400 for the five rungs and 250 for the two references; every
# rung's prior is far from either, so the bound is non-binding throughout.
#
# PREDICTION, recorded in wockner-fit.R before the runs: relative to
# np_wide_total0 the ladder should move cycle_length only ~0.1-0.3 h. A much
# larger move means that reasoning was wrong.
#
# READ THE THREE TOGETHER. A rung that fits as well, samples better, and is
# biologically defensible is the case for fixing b_shape. A rung that buys
# biology at a real predictive cost is a tradeoff to state, not to resolve.
#
# Post-hoc on saved fits, no refitting. Reads seven 55 MB fits one at a time.
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/wockner-bshape-ladder.R 2>&1 | tee _data/bshape-ladder.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({
    library(plasmofit); library(tidyverse); library(rstan); library(loo)
})

LADDER <- tribble(
    ~label,            ~pinned,
    "np_wide_total0",  NA_real_,
    "np_bs50",         50,
    "np_bs84",         84,
    "np_bs100",        100,
    "np_bs150",        150,
    "np_bs250",        250,
    ## max_shape = 1000 from here down, against 400 for the rungs above.
    ## np_bs250_ms1000 is np_bs250 with max_shape ALONE changed, so the gap
    ## between the two is what max_shape is worth on its own -- it sets the
    ## scale of Stan's bounded transform and so can move divergences and step
    ## size without the pin changing at all. Read it before 400 and 600.
    "np_bs250_ms1000", 250,
    "np_bs400",        400,
    "np_bs600",        600,
    "np_wide_both",    NA_real_)

PARS <- c("log10_total0", "R", "cycle_length", "b_shape", "sd_iRBC")
CL_REF <- 45.52   # np_wide_total0's cycle_length, for the age-range column

read_one <- function(cfg) {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    ss <- rstan::summary(f)$summary
    n_post <- sum(sapply(sp, nrow))
    health <- tibble(
        div = sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
        div_pct = 100 * div / n_post,
        max_rhat = max(ss[, "Rhat"], na.rm = TRUE),
        min_ess = min(ss[, "n_eff"], na.rm = TRUE),
        stepsize = mean(sapply(sp, \(x) mean(x[, "stepsize__"]))))
    pars <- map(PARS, \(p) {
        m <- as.matrix(f, pars = p)
        tibble(par = p, mean = mean(colMeans(m)))
    }) |> list_rbind()
    rm(f); invisible(gc())
    list(health = health, pars = pars)
}

cat("=== reading", nrow(LADDER), "fits ===\n")
res <- map(LADDER$label, \(c) { cat("  ", c, "\n"); read_one(c) })
names(res) <- LADDER$label

cat("\n=== 1. sampler health ===\n")
cat("Cells: one row per fit. `div` is divergent transitions after warmup\n")
cat("summed over 4 chains and `div_pct` the same as a percent of the 4000\n")
cat("post-warmup transitions; `max_rhat` the largest split R-hat over all\n")
cat("parameters, `min_ess` the smallest bulk n_eff, `stepsize` the mean\n")
cat("adapted step size over chains. Absolute, not a comparison.\n\n")
map(LADDER$label, \(l) res[[l]]$health |> mutate(label = l, .before = 1)) |>
    list_rbind() |>
    left_join(LADDER, by = "label") |>
    relocate(pinned, .after = label) |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    print(width = Inf)

cat("\n=== 2. posterior means, and what moved ===\n")
cat("Cells: posterior mean of each parameter, averaged over that parameter's\n")
cat("groups (14 grp_init for log10_total0 and b_shape, 13 for R and\n")
cat("cycle_length, 27 for sd_iRBC). `pinned` is where that rung's tight\n")
cat("prior was centred, NA where b_shape was free. `d_cl` is cycle_length\n")
cat("minus np_wide_total0's, in hours -- the quantity the prediction was\n")
cat("about. `age_range_h` is the central 99% of the starting-stage\n")
cat("distribution implied by the fitted b_shape at a 45.5 h cycle.\n\n")
wide <- map(LADDER$label, \(l) res[[l]]$pars |> mutate(label = l)) |>
    list_rbind() |> pivot_wider(names_from = par, values_from = mean) |>
    left_join(LADDER, by = "label") |>
    mutate(d_cl = cycle_length - cycle_length[label == "np_wide_total0"],
           age_range_h = (qbeta(0.995, b_shape, b_shape) -
                          qbeta(0.005, b_shape, b_shape)) * CL_REF)
wide |>
    select(label, pinned, b_shape, age_range_h, cycle_length, d_cl,
           R, log10_total0, sd_iRBC) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(width = Inf)

cat("\n=== 3. predictive cost of pinning ===\n")
cat("Cells: `elpd` is the PSIS-LOO estimate summed over all 1130\n")
cat("observations, higher is better, with its standard error; `k_hi` counts\n")
cat("observations with Pareto k > 0.7. `vs_free` is the PAIRED elpd\n")
cat("difference against np_wide_total0, which leaves b_shape free, and\n")
cat("`vs_bs100` the same against the np_bs100 rung so the rungs can be\n")
cat("compared with each other; both are that row minus the named reference,\n")
cat("summed over the same 1130 observations, POSITIVE MEANS BETTER, with\n")
cat("the standard error of the paired difference beside each. A difference\n")
cat("within about 2 se is not detectable.\n\n")
## The saved LOO object is a list of two groupings (`observation`, `trial`);
## take the observation-level one, as wockner-prior-panel.R's `pick` does.
## Trial-level loo is broken here anyway -- Pareto k > 0.7 for all 13 units.
pick <- function(x) if (inherits(x, "loo")) x else x[["observation"]]
loos <- map(set_names(LADDER$label),
            \(l) pick(read_rds(sprintf("_data/wock-fit-LOO-%s.rds", l))))
## Paired difference computed directly from the pointwise elpd rather than
## through loo_compare(), whose row order depends on the ranking and is easy
## to index wrongly. se_diff is sqrt(n) * sd of the per-observation
## difference, which is what loo_compare reports.
pair <- function(a, b) {
    d <- loos[[a]]$pointwise[, "elpd_loo"] - loos[[b]]$pointwise[, "elpd_loo"]
    c(diff = sum(d), se = sqrt(length(d)) * sd(d))
}
map(LADDER$label, \(l) {
    e <- loos[[l]]$estimates
    v <- if (identical(l, "np_wide_total0")) c(0, 0) else pair(l, "np_wide_total0")
    w <- if (identical(l, "np_bs100")) c(0, 0) else pair(l, "np_bs100")
    tibble(label = l, elpd = e["elpd_loo", "Estimate"], se = e["elpd_loo", "SE"],
           k_hi = sum(loos[[l]]$diagnostics$pareto_k > 0.7),
           vs_free = v[1], se_free = v[2], vs_bs100 = w[1], se_bs100 = w[2])
}) |> list_rbind() |>
    left_join(LADDER, by = "label") |> relocate(pinned, .after = label) |>
    mutate(across(where(is.numeric), \(x) round(x, 2))) |>
    print(width = Inf)

write_rds(list(ladder = LADDER, health = map(res, "health"), means = wide),
          "_data/wock-bshape-ladder.rds")
cat("\nwrote _data/wock-bshape-ladder.rds\n")
