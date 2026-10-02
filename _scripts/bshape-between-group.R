# Do the 14 grp_init groups have different b_shape, or is one value enough?
#
# The ladder pinned all 14 groups at a common centre with a tight prior, so it
# tested FIXED vs ESTIMATED and said nothing about ONE SHARED VALUE vs 14.
# Collapsing `vector[n_grp_init] b_shape` to a scalar is a separate change and
# needs its own evidence.
#
# The test: compare the spread of the 14 posterior means against the typical
# within-group posterior sd. If the groups' estimates are no further apart
# than one group's own uncertainty, the data carry no evidence that they
# differ and a single value loses nothing. Read on the LOG scale, since
# b_shape's prior is lognormal and its posterior is strongly right-skewed.
#
# Reads only the saved summaries. Seconds to run:
#     Rscript --vanilla _scripts/bshape-between-group.R \
#         2>&1 | tee _data/bshape-between-group.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse) })

p <- read_rds("_data/wock-prior-panel.rds")

cat("Cells: one row per fit, over its 14 grp_init groups. `sd_between` is\n")
cat("the sd of the 14 per-group posterior MEANS; `sd_within` is the mean of\n")
cat("the 14 per-group posterior SDs. `ratio` is sd_between / sd_within:\n")
cat("well below 1 means the groups sit closer together than one group's own\n")
cat("uncertainty, so the data do not resolve differences between them.\n")
cat("`range` is the smallest and largest per-group posterior mean. b_shape\n")
cat("is dimensionless; `log` rows are the same quantities computed on\n")
cat("log(b_shape), which is the scale its prior is written on and the scale\n")
cat("a right-skewed positive parameter should be compared on.\n\n")

map(c("no_pool", "np_wide_bshape", "np_wide_both"), \(lab) {
    b <- p$pars |> filter(par == "b_shape", label == lab)
    ## Posterior sd of log(X) approximated by the delta method, sd(X)/mean(X),
    ## because only the natural-scale mean and sd were stored. Adequate here:
    ## it is used to judge an order of magnitude, not to report an interval.
    tibble(fit = lab,
           scale = c("natural", "log"),
           sd_between = c(sd(b$mean), sd(log(b$mean))),
           sd_within = c(mean(b$sd), mean(b$sd / b$mean)),
           lo = c(min(b$mean), log(min(b$mean))),
           hi = c(max(b$mean), log(max(b$mean))))
}) |>
    list_rbind() |>
    mutate(ratio = sd_between / sd_within) |>
    relocate(ratio, .after = sd_within) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf, width = Inf)

cat("\nPer-group posterior means, np_wide_bshape, sorted:\n")
p$pars |> filter(par == "b_shape", label == "np_wide_bshape") |>
    arrange(mean) |> mutate(mean = round(mean, 1), sd = round(sd, 1)) |>
    select(idx, mean, sd) |> print(n = Inf)
