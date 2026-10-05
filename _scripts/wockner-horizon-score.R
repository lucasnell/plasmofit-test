# The horizon ladder: does the hierarchy's advantage come from the TRAINING
# set or from the scoring horizon?
#
# Design A's verdict was mask-dependent -- pooled_cl 3.27 log units behind
# no_pool at the last-third mask, 0.04 behind at the last-quarter. Decomposing
# the third-mask fits put 89% of the deficit on the points the quarter mask
# ALSO holds out, so it is not the extra horizon. What differs is that the
# quarter-mask fits were TRAINED on the 88 deepest points and the third-mask
# fits were not.
#
# THE DESIGN. Three masks, k = 1, 2, 3 points from the end of each series
# (capped at n - 3 so no series keeps fewer than 3), holding out 177, 350 and
# 479 of 1130. They NEST, so the k=1 points are held out by every rung.
# Scoring all three rungs on THAT COMMON WINDOW fixes the scored observations
# and varies only how much was withheld from training -- which is the quantity
# the decomposition identified and which neither Design A mask isolates.
#
# So a trend across rungs here cannot be a horizon effect: the horizon is the
# same 177 points in every row. It can only be the training set.
#
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/wockner-horizon-score.R 2>&1 | tee _data/horizon-score.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

RUNGS <- c(1, 2, 3)
MODELS <- c("no_pool", "pooled_cl")
lme <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }

cfg <- function(k, m) sprintf("daH%d_%s", k, m)

cat("=== reading 6 fits ===\n")
D <- map(set_names(as.character(RUNGS)), \(k)
         read_rds(sprintf("_data/wock-data-%s.rds", cfg(as.integer(k), "no_pool"))))

## every rung must describe the same observations in the same order, or the
## common scoring window is not common
ref <- D[["1"]]
for (k in names(D)) {
    stopifnot(isTRUE(all.equal(D[[k]]$ts, ref$ts)),
              isTRUE(all.equal(D[[k]]$y, ref$y)),
              identical(attr(D[[k]], "levels"), attr(ref, "levels")))
}
score_on <- ref$hold_out == 1L          # the k=1 window, 177 points
## and it must really be held out in every rung
for (k in names(D)) stopifnot(all(D[[k]]$hold_out[score_on] == 1L))
cat("common scoring window:", sum(score_on), "observations, held out in all",
    length(RUNGS), "rungs\n")
cat("held out per rung:", paste(map_int(D, \(d) sum(d$hold_out)), collapse = ", "), "\n")

lpd <- map(set_names(RUNGS), \(k) map(set_names(MODELS), \(m) {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg(k, m)))
    ss <- rstan::summary(f)$summary
    v <- apply(as.matrix(f, pars = "log_lik"), 2, lme)
    h <- list(div = sum(sapply(get_sampler_params(f, inc_warmup = FALSE),
                               \(x) sum(x[, "divergent__"]))),
              max_rhat = max(ss[, "Rhat"], na.rm = TRUE))
    rm(f); invisible(gc()); cat("  ", cfg(k, m), "\n")
    list(lpd = v, health = h)
}))

cat("\n=== 1. sampler health ===\n")
cat("Cells: divergent transitions after warmup summed over 4 chains, and the\n")
cat("largest split R-hat over all parameters. Absolute, not a comparison.\n\n")
map(RUNGS, \(k) map(MODELS, \(m) tibble(rung = k, model = m,
    div = lpd[[as.character(k)]][[m]]$health$div,
    max_rhat = lpd[[as.character(k)]][[m]]$health$max_rhat)) |> list_rbind()) |>
    list_rbind() |> mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    print(width = Inf)

cat("\n=== 2. scored on the COMMON window ===\n")
cat(sprintf("Cells: log pointwise predictive density summed over the SAME %d
observations in every row -- the k=1 points, held out by all three rungs --
computed per observation as log(mean_s exp(log_lik[s,i])), in log units,
higher better. `withheld` is how many observations that rung's fits were
denied in training, out of 1130. `pooled_cl - no_pool` is the paired
difference over those same points with the standard error of the paired
difference; NEGATIVE means collapsing cycle_length predicts worse.\n\n",
            sum(score_on)))

tab <- map(RUNGS, \(k) {
    a <- lpd[[as.character(k)]][["no_pool"]]$lpd[score_on]
    b <- lpd[[as.character(k)]][["pooled_cl"]]$lpd[score_on]
    d <- b - a
    tibble(rung = k, withheld = sum(D[[as.character(k)]]$hold_out),
           no_pool = sum(a), pooled_cl = sum(b),
           diff = sum(d), se = sqrt(length(d)) * sd(d),
           z = sum(d) / (sqrt(length(d)) * sd(d)))
}) |> list_rbind()
tab |> mutate(across(where(is.numeric), \(x) round(x, 2))) |> print(width = Inf)

cat("\n=== verdict ===\n")
cat("  The horizon is identical in every row, so any trend is the training\n")
cat("  set and nothing else.\n")
sl <- coef(lm(diff ~ withheld, data = tab))[2]
cat(sprintf("  slope: %+.4f log units per observation withheld\n", sl))
cat(if (all(abs(tab$z) < 2)) {
    "  FLAT -- pooled_cl is indistinguishable from no_pool at every training\n  size, so Design A's third-mask result does not reproduce here.\n"
} else if (sl < 0 && tab$z[which.max(tab$withheld)] < -2) {
    "  TRAINING-SET EFFECT -- the deficit grows as more is withheld, which is\n  what Design A's mask-dependence implied and is NOT a horizon effect.\n"
} else {
    "  MIXED -- read the rows; the pattern is not a simple monotone trend.\n"
})

write_rds(list(table = tab, n_scored = sum(score_on)),
          "_data/wock-horizon-score.rds")
