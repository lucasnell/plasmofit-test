# Local analysis of wockner-fit-kfold.R's leave-one-trial-out stability check.
#
# Expects wock-fit-kfold-<model>-<trial>.rds to have been copied into _data/.
#
# This is a robustness check, not predictive cross-validation -- see
# wockner-fit-kfold.R's header and claude/CLAUDE.md for why. What matters
# here is whether the population-level cycle-length (and R) estimates hold
# steady regardless of which trial is left out, not which model "wins".


suppressPackageStartupMessages({
    library(tidyverse)
})

res_files <- list.files("_data", "^wock-fit-kfold-.*[.]rds$", full.names = TRUE)
if (length(res_files) == 0) stop("no wock-fit-kfold-*.rds found in _data/")

res <- map(res_files, read_rds)

tbl <- map(res, \(r) {
    tibble(model         = r$model,
           held_trial    = r$held_trial,
           n_grp_cl      = r$n_grp_cl,
           divergences   = r$n_divergent,
           max_rhat      = round(r$max_rhat, 4),
           cl_pop_mean_h = round(r$cl_pop_mean_h, 2),
           cl_lo         = round(unname(r$cl_pop_q[1]), 2),
           cl_hi         = round(unname(r$cl_pop_q[2]), 2),
           R_pop_mean    = round(r$R_pop_mean, 3),
           sigma_cl_mean = round(r$sigma_cl_mean, 4),
           sigma_R_mean  = round(r$sigma_R_mean, 4))
}) |> list_rbind() |> arrange(model, held_trial)

cat("=== per-fold population estimates ===\n")
print(tbl, n = Inf, width = Inf)


# ------------------------------------------------------------------------ #
# Stability across folds: how much does leaving out any one trial move the
# population-level cycle-length estimate? A small spread here means the
# full-data finding isn't riding on any single influential trial.
# ------------------------------------------------------------------------ #

cat("\n=== stability across the 13 folds, by model ===\n")
tbl |>
    summarise(.by = model,
              cl_mean   = round(mean(cl_pop_mean_h), 2),
              cl_spread = round(diff(range(cl_pop_mean_h)), 2),
              cl_sd     = round(sd(cl_pop_mean_h), 3),
              any_divergent = any(divergences > 0),
              max_rhat  = round(max(max_rhat), 4)) |>
    print(width = Inf)

cat("\nsigma_logit_cl across folds (no_pool only; NA for pooled_cl, which has\n",
    "no such parameter):\n", sep = "")
tbl |>
    filter(model == "no_pool") |>
    summarise(sigma_cl_mean_of_folds = round(mean(sigma_cl_mean), 4),
              sigma_cl_max_of_folds  = round(max(sigma_cl_mean), 4)) |>
    print()

cat("\nIf sigma_logit_cl stays small/prior-dominated across all 13 folds (as\n",
    "in the full-data fit) and cl_pop_mean_h is stable across folds for both\n",
    "models, that reinforces the full-data finding. A fold with a much larger\n",
    "sigma_logit_cl, or a cl_pop_mean_h outlier, flags that trial as unusually\n",
    "influential and worth a closer look before trusting the full-data answer.\n",
    sep = "")


# ------------------------------------------------------------------------ #
# no_pool vs pooled_cl agreement, fold by fold. Both models' held-out
# prediction for a *new* trial's cycle length collapses to essentially the
# same population-level quantity (mu_logit_cl for no_pool, the single
# logit_cl for pooled_cl) -- see claude/CLAUDE.md -- so large, systematic
# disagreement between them here would be a surprise worth chasing down
# rather than confirmation that the hierarchy matters.
# ------------------------------------------------------------------------ #

wide <- tbl |>
    select(model, held_trial, cl_pop_mean_h) |>
    pivot_wider(names_from = model, values_from = cl_pop_mean_h)

if (all(c("no_pool", "pooled_cl") %in% names(wide))) {
    cat("\n=== no_pool vs pooled_cl population cycle_length, by held-out trial ===\n")
    wide |>
        mutate(diff_h = round(no_pool - pooled_cl, 3)) |>
        print(n = Inf)
}
