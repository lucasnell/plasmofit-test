# Threads 2 and 3 on the real data, read off one panel of eight fits.
#
# Thread 2 -- is `b_shape`'s prior doing on real data what it does in
# simulation? lognormal(2, 0.5) has median 7.39 against a posterior of 14-19,
# and in simulation widening it removes 0.50 h of cycle-length bias. Two
# contrasts, because the isolated effect and the effect on top of a corrected
# `log10_total0` prior are different questions:
#     no_pool        -> np_wide_bshape   (b_shape alone, old total0 prior)
#     np_wide_total0 -> np_wide_both     (b_shape on top of a fixed total0)
#
# Thread 3 -- every no_pool vs pooled_cl conclusion in this project was drawn
# under the misspecified normal(1, 0.25) prior on `log10_total0`, which costs
# 104.7 elpd. Does the hierarchy conclusion survive a defensible prior? Three
# rungs of the same comparison:
#     no_pool        vs pooled_cl        (as originally run)
#     np_wide_total0 vs pl_wide_total0   (total0 prior corrected)
#     np_wide_both   vs pl_wide_both     (both nuisance priors corrected)
#
# Fits come from `_scripts/wockner-fit.R` tasks 1, 2, 8-13. Runs in a few
# minutes; the fits are ~55 MB each and are read one at a time.
#
#   srun -N 1 -n 1 -c 4 --mem=32G Rscript --vanilla \
#       _scripts/wockner-prior-panel.R 2>&1 | tee _data/prior-panel.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
})

## label -> (config on disk, model, what its priors are). `total0` and
## `bshape` say which nuisance prior is at its default and which is widened,
## which is what both threads are indexed by.
PANEL <- tribble(
    ~label,           ~config,          ~model,      ~total0,  ~bshape,
    "no_pool",        "no_pool",        "no_pool",   "tight",  "tight",
    "np_wide_total0", "np_wide_total0", "no_pool",   "wide",   "tight",
    "np_wide_bshape", "np_wide_bshape", "no_pool",   "tight",  "wide",
    "np_wide_both",   "np_wide_both",   "no_pool",   "wide",   "wide",
    "np_anchor",      "np_anchor",      "no_pool",   "anchor", "tight",
    "pooled_cl",      "pooled_cl",      "pooled_cl", "tight",  "tight",
    "pl_wide_total0", "pl_wide_total0", "pooled_cl", "wide",   "tight",
    "pl_wide_both",   "pl_wide_both",   "pooled_cl", "wide",   "wide")

PARS <- c("log10_total0", "R", "cycle_length", "b_shape", "sd_iRBC")

fit_path <- \(cfg) sprintf("_data/wock-fit-%s.rds", cfg)
loo_path <- \(cfg) sprintf("_data/wock-fit-LOO-%s.rds", cfg)

have <- file.exists(fit_path(PANEL$config))
if (!all(have)) {
    cat("MISSING fits, these rows are dropped:\n")
    cat(sprintf("  %s\n", PANEL$label[!have]), sep = "")
    cat("\n")
}
PANEL <- PANEL[have, ]
stopifnot(nrow(PANEL) >= 2)


# ------------------------------------------------------------------------ #
# One pass over the fits: sampler health and per-group posterior summaries.
# Kept in memory as small tibbles so the 55 MB fit objects can be dropped.
# ------------------------------------------------------------------------ #

read_one <- function(cfg) {
    f <- read_rds(fit_path(cfg))
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    ss <- rstan::summary(f)$summary
    n_post <- sum(sapply(sp, nrow))
    health <- tibble(
        div = sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
        div_pct = 100 * div / n_post,
        treedepth = sum(sapply(sp, \(x) sum(x[, "treedepth__"] >= 10))),
        max_rhat = max(ss[, "Rhat"], na.rm = TRUE),
        min_ess = min(ss[, "n_eff"], na.rm = TRUE),
        stepsize = mean(sapply(sp, \(x) mean(x[, "stepsize__"]))))
    pars <- map(PARS, \(p) {
        m <- as.matrix(f, pars = p)
        tibble(par = p, idx = seq_len(ncol(m)),
               mean = colMeans(m), sd = apply(m, 2, sd),
               ## MCSE of the posterior mean, for judging whether a shift
               ## between two fits is bigger than Monte Carlo noise.
               mcse = apply(m, 2, sd) / sqrt(ss[colnames(m), "n_eff"]))
    }) |> list_rbind()
    rm(f); invisible(gc())
    list(health = health, pars = pars)
}

cat("=== reading", nrow(PANEL), "fits ===\n")
res <- map(PANEL$config, \(cfg) { cat("  ", cfg, "\n"); read_one(cfg) })
names(res) <- PANEL$label

health <- map2(res, PANEL$label, \(r, lab) mutate(r$health, label = lab)) |>
    list_rbind() |> relocate(label)
pars <- map2(res, PANEL$label, \(r, lab) mutate(r$pars, label = lab)) |>
    list_rbind() |> relocate(label)


# ------------------------------------------------------------------------ #
# 1. Sampler health. Everything below is void for any row that fails this,
#    so it is printed first and gated on explicitly.
# ------------------------------------------------------------------------ #

cat("\n=== 1. sampler health ===\n")
cat("Cells: one row per fit. `div` is the count of divergent transitions\n")
cat("after warmup summed over 4 chains (`div_pct` as a percent of the 4000\n")
cat("post-warmup transitions), `treedepth` the count hitting max_treedepth,\n")
cat("`max_rhat` the largest split R-hat over all parameters, `min_ess` the\n")
cat("smallest bulk n_eff over all parameters, `stepsize` the mean adapted\n")
cat("step size over chains. Not a comparison -- no baseline.\n\n")
health |> mutate(across(where(is.numeric), \(x) signif(x, 4))) |>
    print(width = Inf)

ok <- health$max_rhat < 1.05
if (!all(ok)) {
    cat("\n  FAILED CONVERGENCE, dropped from everything below:\n")
    cat(sprintf("    %s (max R-hat %.4f)\n",
                health$label[!ok], health$max_rhat[!ok]), sep = "")
    keep <- health$label[ok]
    PANEL <- PANEL[PANEL$label %in% keep, ]
    pars <- pars[pars$label %in% keep, ]
    health <- health[ok, ]
} else {
    cat("\n  all max R-hat < 1.05\n")
}
avail <- PANEL$label


# ------------------------------------------------------------------------ #
# 2. Parameters, averaged over groups.
# ------------------------------------------------------------------------ #

wide_mean <- pars |>
    summarise(mean = mean(mean), .by = c(label, par)) |>
    pivot_wider(names_from = par, values_from = mean) |>
    left_join(PANEL |> select(label, model, total0, bshape), by = "label") |>
    relocate(label, model, total0, bshape) |>
    arrange(match(label, PANEL$label))

cat("\n=== 2. posterior means, averaged over groups ===\n")
cat("Cells: the posterior mean of each parameter, then averaged over that\n")
cat("parameter's groups within a fit -- 14 grp_init for log10_total0, 13\n")
cat("trials for R and cycle_length, 1 for b_shape, 27 for sd_iRBC. Units:\n")
cat("log10_total0 in log10 iRBC per ml, R and b_shape dimensionless,\n")
cat("cycle_length in hours, sd_iRBC in log10 units. Not a comparison --\n")
cat("each row stands alone; `total0` and `bshape` say which nuisance prior\n")
cat("is at its default (\"tight\"), widened (\"wide\"), or replaced by the\n")
cat("inoculum anchor.\n\n")
wide_mean |> mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(width = Inf)


# ------------------------------------------------------------------------ #
# 3. Predictive comparison.
# ------------------------------------------------------------------------ #

pick <- function(x) if (inherits(x, "loo")) x else x[[1]]
has_loo <- file.exists(loo_path(PANEL$config))
loos <- set_names(map(PANEL$config[has_loo], \(c) pick(read_rds(loo_path(c)))),
                  PANEL$label[has_loo])

cat("\n=== 3. predictive comparison (loo, per observation) ===\n")
if (!all(has_loo))
    cat("  no loo for:", paste(PANEL$label[!has_loo], collapse = ", "), "\n")
cat("Cells: `elpd` is the PSIS-LOO estimate of expected log pointwise\n")
cat("predictive density summed over all 1130 observations, higher is\n")
cat("better, with its standard error; `k_hi` counts observations with\n")
cat("Pareto k > 0.7, where the importance-sampling approximation is\n")
cat("unreliable. Absolute, not a comparison.\n\n")
elpd <- imap(loos, \(l, nm) tibble(
    label = nm,
    elpd = l$estimates["elpd_loo", 1], se = l$estimates["elpd_loo", 2],
    p_loo = l$estimates["p_loo", 1],
    k_hi = sum(l$diagnostics$pareto_k > 0.7))) |> list_rbind()
elpd |> mutate(across(where(is.numeric), \(x) round(x, 1))) |> print(n = Inf)

if (length(loos) >= 2) {
    cat("\nCells: `elpd_diff` is each fit's elpd minus the best fit's elpd\n")
    cat("(so 0 at the top, negative below, in elpd units summed over 1130\n")
    cat("observations), and `se_diff` the standard error of that paired\n")
    cat("difference -- paired, since it is computed per observation on the\n")
    cat("same data. A difference of a few se is a real predictive gain.\n\n")
    cmpl <- loo::loo_compare(loos)
    cp <- as.data.frame(cmpl)
    num <- vapply(cp, is.numeric, logical(1))
    cp[num] <- lapply(cp[num], \(x) round(x, 2))
    print(cp[, intersect(c("model", "elpd_diff", "se_diff"), names(cp)),
             drop = FALSE])
} else cmpl <- NULL

## Paired elpd difference between two named fits, with its own se. This is
## the quantity each thread asks for; loo_compare only ranks against the best.
elpd_pair <- function(a, b) {
    if (!all(c(a, b) %in% names(loos))) return(tibble(from = a, to = b))
    pa <- loos[[a]]$pointwise[, "elpd_loo"]
    pb <- loos[[b]]$pointwise[, "elpd_loo"]
    d <- pb - pa
    tibble(from = a, to = b, elpd_diff = sum(d),
           se_diff = sqrt(length(d)) * sd(d),
           z = sum(d) / (sqrt(length(d)) * sd(d)))
}


# ------------------------------------------------------------------------ #
# 4. Thread 2: the b_shape prior.
# ------------------------------------------------------------------------ #

## Paired shift in a parameter between two fits, over that parameter's own
## groups. Paired because both fits see the same data and the same groups in
## the same order, so group i in one matches group i in the other.
shift <- function(a, b, parname) {
    x <- pars[pars$label %in% c(a, b) & pars$par == parname, ]
    if (n_distinct(x$label) < 2) return(tibble(par = parname, from = a, to = b))
    w <- x |> select(label, idx, mean) |>
        pivot_wider(names_from = label, values_from = mean)
    d <- w[[b]] - w[[a]]
    ## MCSE of the paired difference: the two fits are independent, so their
    ## Monte Carlo errors add in quadrature.
    m <- x |> select(label, idx, mcse) |>
        pivot_wider(names_from = label, values_from = mcse)
    se <- sqrt(m[[a]]^2 + m[[b]]^2)
    tibble(par = parname, from = a, to = b, n_grp = length(d),
           mean_from = mean(w[[a]]), mean_to = mean(w[[b]]),
           mean_shift = mean(d), max_abs_shift = max(abs(d)),
           mean_z = mean(d / se))
}

T2 <- list(c("no_pool", "np_wide_bshape"),
           c("np_wide_total0", "np_wide_both"))

cat("\n\n=== 4. THREAD 2: the b_shape prior on real data ===\n")
cat("Two contrasts. `no_pool -> np_wide_bshape` widens sd_log_b_shape from\n")
cat("0.5 to 1.5 with everything else at its default. `np_wide_total0 ->\n")
cat("np_wide_both` does the same on top of the corrected log10_total0\n")
cat("prior, so it asks whether the b_shape effect survives once the prior\n")
cat("known to be misspecified is fixed.\n")
cat("\nCells: `mean_from` and `mean_to` are the posterior mean averaged over\n")
cat("the parameter's groups in each fit; `mean_shift` is `to` minus\n")
cat("`from`, averaged over groups, paired by group; `max_abs_shift` the\n")
cat("largest absolute per-group shift; `mean_z` the paired shift divided by\n")
cat("its Monte Carlo standard error, averaged over groups, so |z| of order\n")
cat("1 is indistinguishable from sampling noise. Units are each\n")
cat("parameter's own -- hours for cycle_length, log10 units for\n")
cat("log10_total0 and sd_iRBC, dimensionless for R and b_shape.\n\n")
t2 <- map(T2, \(p) map(PARS, \(q) shift(p[1], p[2], q)) |> list_rbind()) |>
    list_rbind()
t2 |> mutate(across(where(is.numeric), \(x) round(x, 4))) |> print(width = Inf, n = Inf)

cat("\nCells: paired PSIS-LOO difference, `to` minus `from`, summed over the\n")
cat("1130 observations both fits share; positive favours `to`. `z` is that\n")
cat("difference over its standard error.\n\n")
map(T2, \(p) elpd_pair(p[1], p[2])) |> list_rbind() |>
    mutate(across(where(is.numeric), \(x) round(x, 2))) |> print(width = Inf)

cat("\nRead: in simulation, widening b_shape's prior removed 0.501 h of\n")
cat("cycle-length bias (-0.531, -0.452, -0.519 over three paired\n")
cat("replicates) and widening log10_total0's removed none. On real data\n")
cat("there is no truth to compare against, so the cycle_length row above is\n")
cat("a shift, not a bias: it says how much of the reported cycle length is\n")
cat("the prior's doing, which bounds how much of it can be biological.\n")


# ------------------------------------------------------------------------ #
# 5. Thread 3: does the hierarchy conclusion survive a defensible prior?
# ------------------------------------------------------------------------ #

T3 <- list(c("no_pool", "pooled_cl"),
           c("np_wide_total0", "pl_wide_total0"),
           c("np_wide_both", "pl_wide_both"))

cat("\n\n=== 5. THREAD 3: no_pool vs pooled_cl at three priors ===\n")
cat("Same comparison run three times, at successively more defensible\n")
cat("nuisance priors: as originally run, with log10_total0 corrected, and\n")
cat("with both corrected. If the original conclusion is an artefact of the\n")
cat("misspecified prior, the sign or size of these rows will move.\n")
cat("\nCells: as in section 4, but `from` is always the no_pool member of\n")
cat("the pair and `to` the pooled_cl member, so `mean_shift` is pooled_cl\n")
cat("minus no_pool. cycle_length is the one that matters: it is the\n")
cat("pooling offset, in hours, averaged over 13 trials and paired by\n")
cat("trial.\n\n")
t3 <- map(T3, \(p) map(PARS, \(q) shift(p[1], p[2], q)) |> list_rbind()) |>
    list_rbind()
t3 |> mutate(across(where(is.numeric), \(x) round(x, 4))) |> print(width = Inf, n = Inf)

cat("\nCells: paired PSIS-LOO difference, pooled_cl minus no_pool, summed\n")
cat("over 1130 observations; positive favours pooled_cl.\n\n")
map(T3, \(p) elpd_pair(p[1], p[2])) |> list_rbind() |>
    mutate(across(where(is.numeric), \(x) round(x, 2))) |> print(width = Inf)

## The pooling offset on its own, across the three prior settings, since it
## is the single number the hierarchy conclusion rests on.
cat("\nCells: the cycle-length pooling offset only -- posterior mean\n")
cat("cycle_length under pooled_cl minus under no_pool, in hours, averaged\n")
cat("over 13 trials, paired by trial. One row per prior setting.\n\n")
t3 |> filter(par == "cycle_length") |>
    transmute(priors = case_when(from == "no_pool" ~ "as originally run",
                                 from == "np_wide_total0" ~ "log10_total0 corrected",
                                 TRUE ~ "both corrected"),
              no_pool = round(mean_from, 3), pooled_cl = round(mean_to, 3),
              offset_h = round(mean_shift, 3), mean_z = round(mean_z, 2)) |>
    print(width = Inf)


out <- list(panel = PANEL, health = health, pars = pars,
            means = wide_mean, elpd = elpd, loo_compare = cmpl,
            thread2 = t2, thread3 = t3,
            thread2_elpd = map(T2, \(p) elpd_pair(p[1], p[2])) |> list_rbind(),
            thread3_elpd = map(T3, \(p) elpd_pair(p[1], p[2])) |> list_rbind())
write_rds(out, "_data/wock-prior-panel.rds")
cat("\nwrote _data/wock-prior-panel.rds\n")
