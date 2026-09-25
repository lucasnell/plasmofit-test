# Does the panel's `b_shape` posterior mean mean anything?
#
# Section 2 of `_scripts/wockner-prior-panel.R` reports `b_shape` moving from
# a posterior mean of ~14.9 to ~65 when `sd_log_b_shape` goes 0.5 -> 1.5. A
# mean is a bad summary of a weakly identified positive parameter, so this
# reads the posterior sd already stored in `_data/wock-prior-panel.rds` and
# asks whether the wide-prior mean sits on a posterior that is merely
# relocated or one that is spread out to the width of its own prior.
#
# Reads only the 87 KB summary object, not the 55 MB fits. Seconds to run:
#     Rscript --vanilla _scripts/panel-bshape-spread.R \
#         2>&1 | tee _data/panel-bshape-spread.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse) })

p <- read_rds("_data/wock-prior-panel.rds")

## Prior coefficient of variation of a lognormal(mu, sigma) is
## sqrt(exp(sigma^2) - 1), independent of mu. The two settings in the panel:
prior_cv <- \(s) sqrt(exp(s^2) - 1)
cat(sprintf("prior CV of b_shape: tight (sd_log = 0.5) = %.3f,  wide (1.5) = %.3f\n\n",
            prior_cv(0.5), prior_cv(1.5)))

cat("Cells: one row per fit. `mean` and `sd` are the posterior mean and\n")
cat("posterior sd of b_shape, averaged over its 14 grp_init groups; `CV` is\n")
cat("the mean over groups of each group's own sd/mean, so it is a posterior\n")
cat("width relative to location and is comparable to the prior CV above.\n")
cat("`bshape` says whether that fit's b_shape prior is at its default\n")
cat("(tight) or widened. b_shape is dimensionless.\n\n")

p$pars |>
    filter(par == "b_shape") |>
    group_by(label) |>
    summarise(mean = mean(mean), sd = mean(sd), CV = mean(sd / mean),
              .groups = "drop") |>
    left_join(select(p$panel, label, bshape), by = "label") |>
    arrange(bshape, label) |>
    mutate(across(c(mean, sd, CV), \(x) round(x, 3))) |>
    print(n = Inf)

cat("\nCells: as above, for cycle_length over its 13 grp_cl groups, in hours.\n")
cat("`sd` is the within-fit posterior sd, i.e. how well one trial's cycle\n")
cat("length is pinned down, not the spread across trials.\n\n")

p$pars |>
    filter(par == "cycle_length") |>
    group_by(label) |>
    summarise(mean = mean(mean), sd = mean(sd), .groups = "drop") |>
    left_join(select(p$panel, label, total0, bshape), by = "label") |>
    arrange(total0, bshape) |>
    mutate(across(c(mean, sd), \(x) round(x, 3))) |>
    print(n = Inf)
