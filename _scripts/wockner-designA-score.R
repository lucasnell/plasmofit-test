# Design A: score the four model variants on the window none of them saw.
#
# Every no_pool vs pooled_cl verdict in this project has rested on a
# comparison its own documentation calls weak -- observation-level loo, where
# a held-out point is pinned by its neighbours whatever the model -- or on
# trial-level loo, which is broken here (Pareto k > 0.7 for all 13 units, in
# every model). Design A masks the last third of every series, fits, and
# scores the held-out window directly. No importance sampling is involved and
# none of its diagnostics apply: these points were genuinely not fitted.
#
# THE QUANTITY. For each held-out observation i, the log pointwise predictive
# density is lpd_i = log( mean_s exp(log_lik[s, i]) ) over posterior draws s,
# computed by log-sum-exp. Held-out elpd is the sum of lpd_i over the held-out
# points only. generated quantities computes log_lik for every observation,
# fitted or not, so the in-sample sum is available too and is reported beside
# it -- but the held-out column is the comparison. Higher is better.
#
# Differences are PAIRED: every model is scored on the same observations, so
# the standard error of a difference is sqrt(n) * sd of the per-observation
# difference, which is much smaller than the se of either total.
#
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/wockner-designA-score.R 2>&1 | tee _data/designA-score.log
#
# PREFIX=daQ_ scores the quarter-mask arm instead of the third.

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

PREFIX <- Sys.getenv("PREFIX", "daA_")
MODELS <- c("no_pool", "pooled_cl", "pooled_R", "pooled_both")
REF <- "no_pool"

lme <- function(x) {        # log mean exp, stable
    m <- max(x)
    m + log(mean(exp(x - m)))
}

cat("=== reading", length(MODELS), "fits, prefix", PREFIX, "===\n")
res <- map(set_names(MODELS), \(m) {
    cfg <- paste0(PREFIX, m)
    d <- read_rds(sprintf("_data/wock-data-%s.rds", cfg))
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    ll <- as.matrix(f, pars = "log_lik")
    sp <- get_sampler_params(f, inc_warmup = FALSE)
    ss <- rstan::summary(f)$summary
    health <- tibble(div = sum(sapply(sp, \(x) sum(x[, "divergent__"]))),
                     max_rhat = max(ss[, "Rhat"], na.rm = TRUE),
                     min_ess = min(ss[, "n_eff"], na.rm = TRUE))
    rm(f); invisible(gc())
    stopifnot(ncol(ll) == d$n_total_obs, length(d$hold_out) == d$n_total_obs)
    cat(sprintf("  %-22s draws %d, obs %d, held out %d\n", cfg, nrow(ll),
                ncol(ll), sum(d$hold_out)))
    list(lpd = apply(ll, 2, lme), ho = d$hold_out == 1L, health = health)
})

ho <- res[[1]]$ho
stopifnot(all(map_lgl(res, \(r) identical(r$ho, ho))))   # same mask everywhere

cat("\n=== 1. sampler health ===\n")
cat("Cells: divergent transitions after warmup summed over 4 chains, largest\n")
cat("split R-hat, smallest bulk n_eff, over all parameters. Absolute.\n\n")
map(MODELS, \(m) res[[m]]$health |> mutate(model = m, .before = 1)) |>
    list_rbind() |> mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    print(width = Inf)

cat("\n=== 2. held-out elpd ===\n")
cat(sprintf("Cells: `held_out` is the log pointwise predictive density summed over
the %d observations NONE of these models were fitted on; `in_sample` the
same over the %d that were. Higher is better, both in log units.
`vs_%s` is the PAIRED held-out difference against %s, that model
minus the reference over the same held-out points, with the standard
error of the paired difference. Positive favours that model.\n\n",
            sum(ho), sum(!ho), REF, REF))

tab <- map(MODELS, \(m) {
    d <- res[[m]]$lpd[ho] - res[[REF]]$lpd[ho]
    tibble(model = m,
           held_out = sum(res[[m]]$lpd[ho]),
           in_sample = sum(res[[m]]$lpd[!ho]),
           vs_ref = sum(d),
           se_diff = if (m == REF) 0 else sqrt(length(d)) * sd(d),
           z = if (m == REF) NA_real_ else sum(d) / (sqrt(length(d)) * sd(d)))
}) |> list_rbind()
tab |> mutate(across(where(is.numeric), \(x) round(x, 2))) |> print(width = Inf)

cat("\n=== convergence gate ===\n")
bad <- map_chr(MODELS, \(m) if (res[[m]]$health$max_rhat >= 1.05) m else NA_character_)
bad <- bad[!is.na(bad)]
if (length(bad)) {
    cat("  VOID (max R-hat >= 1.05):", paste(bad, collapse = ", "),
        "-- everything about these rows is void, per claude/gotchas.md\n")
} else cat("  all models clear R-hat < 1.05\n")

cat("\n=== verdict ===\n")
## NOTE ON DIRECTION, because it is easy to state backwards: `no_pool` is the
## model that KEEPS the hierarchy -- it pools nothing -- and `pooled_cl` is
## the one that collapses cycle_length to a single value. So no_pool scoring
## HIGHER means the hierarchy earns its keep.
hier <- tab |> filter(model == "pooled_cl")
cat(sprintf("  pooled_cl (hierarchy removed) vs no_pool (hierarchy kept),\n  held-out: %+.2f (se %.2f, z %+.1f)\n",
            hier$vs_ref, hier$se_diff, hier$z))
cat(if (abs(hier$z) < 2) {
    "  INDISTINGUISHABLE -- the hierarchy neither earns its keep nor costs,\n  now on a test that is not broken.\n"
} else if (hier$vs_ref < 0) {
    "  THE HIERARCHY EARNS ITS KEEP -- collapsing cycle_length to one value\n  costs held-out predictive accuracy. This REVERSES the lean from the\n  weak comparisons, which put the two indistinguishable.\n"
} else {
    "  POOLING WINS -- collapsing cycle_length improves held-out prediction.\n"
})
if (length(bad)) cat("  (read none of the above for the voided models)\n")

write_rds(list(prefix = PREFIX, table = tab, n_held_out = sum(ho)),
          sprintf("_data/wock-designA-score-%s.rds", sub("_$", "", PREFIX)))
cat("\nwrote _data/wock-designA-score-", sub("_$", "", PREFIX), ".rds\n", sep = "")
