# Why do the third- and quarter-mask verdicts disagree?
#
# Design A at the last THIRD put pooled_cl 3.27 log units behind no_pool on
# held-out data (z -4.5). At the last QUARTER the same comparison is -0.04
# (z -0.1). Both fractions were fixed before any fit ran, so neither can be
# preferred for reading better; the disagreement is the result and it needs
# an explanation rather than a choice.
#
# The two masks NEST: both hold out the final k observations of each series,
# with k = floor(n/3) and floor(n/4), at least 1. So the quarter's 218 points
# are a subset of the third's 306, and the 88 that differ are the ones held
# out ONLY by the third -- which, because both take from the end, are the
# points deeper into each series' tail. Verified, not assumed, in the
# subset check below.
#
# That makes the third-mask fits decomposable WITHOUT any refitting: score
# them separately on the 218 shared points and on the 88 deeper ones. If the
# hierarchy's advantage is concentrated in the deeper points, the mechanism
# is the one Design A was built around -- a cycle-length error accumulating
# as phase drift, which needs horizon to show itself -- and the quarter mask
# simply does not reach far enough. If the advantage is spread evenly, that
# explanation fails and the disagreement is something else.
#
# Within-fit decomposition only. The quarter-mask FITS saw more data (912
# points against 824) and are not comparable to the third-mask fits
# observation by observation; they are reported alongside, not differenced.
#
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/designA-horizon.R 2>&1 | tee _data/designA-horizon.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

MODELS <- c("no_pool", "pooled_cl", "pooled_R", "pooled_both")
REF <- "no_pool"
lme <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }

a <- read_rds("_data/wock-data-daA_no_pool.rds")
q <- read_rds("_data/wock-data-daQ_no_pool.rds")
stopifnot(isTRUE(all.equal(a$ts, q$ts)), isTRUE(all.equal(a$y, q$y)),
          identical(attr(a, "levels"), attr(q, "levels")),
          all(q$hold_out[a$hold_out == 0] == 0))     # quarter nests in third

shared <- a$hold_out == 1 & q$hold_out == 1    # 218, held out by both
deeper <- a$hold_out == 1 & q$hold_out == 0    # 88, held out only by the third

cat("=== reading the third-mask fits ===\n")
lpd <- map(set_names(MODELS), \(m) {
    f <- read_rds(sprintf("_data/wock-fit-daA_%s.rds", m))
    v <- apply(as.matrix(f, pars = "log_lik"), 2, lme)
    rm(f); invisible(gc()); cat("  ", m, "\n"); v
})

cat("\nCells: log pointwise predictive density summed over the stated subset
of held-out observations, in log units, higher better. All four models are
the THIRD-mask fits, so every column scores the same posteriors on
different subsets of the window they did not see. `vs no_pool` is the
paired difference over that subset with the standard error of the paired
difference; negative means the row predicts worse than no_pool.
`shared` is the 218 points the quarter mask also holds out; `deeper` the
88 held out only by the third, which lie further into each series' tail.\n\n")

tab <- map(MODELS, \(m) {
    ds <- lpd[[m]][shared] - lpd[[REF]][shared]
    dd <- lpd[[m]][deeper] - lpd[[REF]][deeper]
    tibble(model = m,
           shared_n = sum(shared), shared_diff = sum(ds),
           shared_se = if (m == REF) 0 else sqrt(length(ds)) * sd(ds),
           deeper_n = sum(deeper), deeper_diff = sum(dd),
           deeper_se = if (m == REF) 0 else sqrt(length(dd)) * sd(dd),
           total_diff = sum(ds) + sum(dd))
}) |> list_rbind()
tab |> mutate(across(where(is.numeric), \(x) round(x, 2))) |> print(width = Inf)

cat("\n=== per-observation rate ===\n")
cat("Cells: the same paired differences divided by the number of points in
each subset, so the two are comparable despite 218 against 88. Log units
per held-out observation; negative means worse than no_pool.\n\n")
tab |> transmute(model,
                 per_obs_shared = shared_diff / shared_n,
                 per_obs_deeper = deeper_diff / deeper_n,
                 ratio = (deeper_diff / deeper_n) / (shared_diff / shared_n)) |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |> print(width = Inf)

h <- tab |> filter(model == "pooled_cl")
cat(sprintf("\npooled_cl: %.2f of its %.2f total comes from the %d deeper points (%.0f%%),\n",
            h$deeper_diff, h$total_diff, h$deeper_n,
            100 * h$deeper_diff / h$total_diff))
cat(sprintf("which are %.0f%% of the held-out set.\n",
            100 * h$deeper_n / (h$shared_n + h$deeper_n)))
