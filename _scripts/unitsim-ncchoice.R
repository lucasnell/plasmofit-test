## What do the simulations say about CHOOSING a fit n_c from outside the data?
## (2026-11-09, Lucas: "check against simulations first", before any decision
## to refit the real data at n_c 512 or 1024 on Park et al.'s figure of
## 450-1009.)
##
## The question is asymmetric and the existing fits answer only half of it:
##   A. fit n_c ABOVE the truth  — gen 384 fitted at 288 / 384 / 512, the
##      fixed-pin grid (truth-centred cycle_length prior, b_shape fixed).
##   B. fit n_c BELOW the truth  — gen 1024 fitted at 384, the main study and
##      its arms (production prior, b_shape estimated or fixed at the truth).
## There is NO gen 1024 / fit 512 cell anywhere, so the combination Park's
## figure would actually put us in is not simulated. Stated, not papered over.
##
## Cells: bias = posterior mean cycle_length minus the truth, hours, positive
## = long, mean over passing fits of 60 per cell. Paired = per replicate on
## identical data, both fits passing: mean (se) over pairs.
##
##   Rscript --vanilla _scripts/unitsim-ncchoice.R > _data/unitsim-read-ncchoice-<date>.txt

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 150)
D <- "_data/unitsim/study/"
x <- bind_rows(lapply(Sys.glob(paste0(D, "sum-*.rds")), readRDS))
for (v in c("prior", "truth_set", "pin")) if (!v %in% names(x)) x[[v]] <- NA_character_
if (!"fit_nc" %in% names(x)) x$fit_nc <- NA_integer_
x <- x |> mutate(prior = coalesce(prior, "prod"), truth_set = coalesce(truth_set, "main"),
                 pin = coalesce(pin, "truth"), fit_nc = coalesce(fit_nc, 384L),
                 key = paste0("w", window))
cat(nrow(x), "merged fits;", sum(x$pass), "pass.\n\n")

pair <- function(a, b, by) {
    j <- inner_join(a, b, by = by, suffix = c("_a", "_b"))
    j |> summarise(pairs = n(), mean = mean(bias_a - bias_b),
                   se = sd(bias_a - bias_b) / sqrt(n()),
                   negative = sum(bias_a - bias_b < 0), .groups = "drop")
}

cat("A. FIT n_c ABOVE THE TRUTH. Truth gen n_c 384, b_shape fixed at the truth\n",
    "   (20.10), cycle_length prior centred at the truth; the fixed-pin grid.\n",
    "   Moving the fit from 384 to 512 is the move Park's figure would justify.\n\n", sep = "")
A <- filter(x, truth_set == "main", gen_nc == 384L, method == "fix",
            prior == "truth", pin == "truth", pass)
A |> group_by(key, fit_nc) |>
    summarise(n_pass = n(), mean_bias = round(mean(bias), 3),
              se = round(sd(bias) / sqrt(n()), 3), cover80 = round(mean(in_80), 3),
              .groups = "drop") |> as.data.frame() |> print(row.names = FALSE)
cat("\n   Paired, identical data: fit 512 minus fit 384, and fit 288 minus fit 384.\n\n")
for (nc in c(512L, 288L)) {
    p <- pair(filter(A, fit_nc == nc), filter(A, fit_nc == 384L), c("key", "rep")) |>
        mutate(contrast = sprintf("fit %d - fit 384", nc))
    print(as.data.frame(mutate(p, across(where(is.double), \(v) round(v, 3)))), row.names = FALSE)
}

cat("\nB. FIT n_c BELOW THE TRUTH. Gen n_c 1024 fitted at 384, main truth,\n",
    "   production prior; b_shape estimated or fixed at the truth. 'mapped' is\n",
    "   the main study (b_shape lowered to hold the total spread at 168 h);\n",
    "   'not mapped' is the bconst arm (b_shape kept at 20.10). Gen 384 at the\n",
    "   same settings is shown as the correctly specified reference.\n\n", sep = "")
B <- filter(x, truth_set %in% c("main", "bconst"), prior == "prod",
            fit_nc == 384L, pass) |>
    mutate(arm = case_when(gen_nc == 384L ~ "gen 384 (correct)",
                           truth_set == "bconst" ~ "gen 1024, b_shape not mapped",
                           TRUE ~ "gen 1024, b_shape mapped"))
B |> group_by(arm, key, method) |>
    summarise(n_pass = n(), mean_bias = round(mean(bias), 3),
              se = round(sd(bias) / sqrt(n()), 3), cover80 = round(mean(in_80), 3),
              .groups = "drop") |> as.data.frame() |> print(row.names = FALSE)
cat("\n   Paired, same noise draws: gen 1024 minus gen 384, both fitted at 384.\n\n")
for (ts in c("main", "bconst")) {
    hi <- filter(B, gen_nc == 1024L, truth_set == ts)
    if (nrow(hi) == 0) next
    p <- pair(hi, filter(B, gen_nc == 384L), c("key", "method", "rep")) |>
        mutate(arm = if (ts == "bconst") "not mapped" else "mapped")
    print(as.data.frame(mutate(p, across(where(is.double), \(v) round(v, 3)))), row.names = FALSE)
}

cat("\nC. THE ASYMMETRY, which is the decision. Mean |bias| in hours over the\n",
    "   cells above, as the cost of each kind of n_c error.\n\n", sep = "")
tibble(
    situation = c("fit 512, truth 384 (overshoot by 1.33x)",
                  "fit 384, truth 384 (correct)",
                  "fit 288, truth 384 (undershoot by 0.75x)",
                  "fit 384, truth 1024 mapped (undershoot by 0.375x)",
                  "fit 384, truth 1024 not mapped (undershoot by 0.375x)"),
    mean_abs_bias = c(
        mean(abs(filter(A, fit_nc == 512L)$bias)),
        mean(abs(filter(A, fit_nc == 384L)$bias)),
        mean(abs(filter(A, fit_nc == 288L)$bias)),
        mean(abs(filter(B, gen_nc == 1024L, truth_set == "main")$bias)),
        mean(abs(filter(B, gen_nc == 1024L, truth_set == "bconst")$bias)))) |>
    mutate(mean_abs_bias = round(mean_abs_bias, 3)) |> as.data.frame() |> print(row.names = FALSE)
cat("\n   NOTE the priors differ between A (centred at the truth) and B\n",
    "   (production, centred 48 h), so A and B are NOT directly comparable in\n",
    "   level; the production-prior residual at gen 384 / fit 384 is in B's own\n",
    "   first rows. The comparison that stands is WITHIN each block.\n", sep = "")
