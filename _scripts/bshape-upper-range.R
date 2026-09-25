# Is fixing b_shape near 100 consistent with these data, and how much would it
# change what the model can represent?
#
# b_shape is the shape of a SYMMETRIC Beta(b_shape, b_shape) over position in
# the cycle (inst/stan/functions/plasmofit.stan, beta_starts), discretised
# into n_c = 96 stage cells and declared `<lower=2, upper=max_shape>` with
# max_shape = 250. Two things follow that the posterior mean alone does not
# show:
#
#   1. whether the widened-prior posterior is running into the 250 bound, in
#      which case its mean of 65.4 is partly an artefact of where the model
#      was cut off rather than a statement about the data;
#   2. how much observable difference there is between b_shape = 65, 100 and
#      250, since the width of the starting-stage distribution falls off as
#      b_shape^(-1/2) and the claim "changing b_shape does little at high
#      values" needs a number.
#
# For a symmetric Beta(a, a) on [0, 1] the sd is 1 / (2 * sqrt(2a + 1)), so
# the starting-stage spread in hours is cycle_length / (2 * sqrt(2a + 1)).
#
# Post-hoc on saved fits. Reads three 55 MB fits one at a time.
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/bshape-upper-range.R 2>&1 | tee _data/bshape-upper-range.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

MAX_SHAPE <- 250
CL <- 45.3   # population cycle_length, h, for converting spread to hours

cat("=== 1. what b_shape means in hours ===\n\n")
cat("Cells: `stage_sd_h` is the sd of the starting-stage distribution in\n")
cat("hours, cycle_length / (2 * sqrt(2 * b_shape + 1)) at cycle_length =\n")
cat(sprintf("%.1f h; `cells` is that in units of one stage compartment\n", CL))
cat("(n_c = 96, so 0.472 h each); `d_sd_per_10` is how much stage_sd_h\n")
cat("changes for a further +10 in b_shape, in hours. Deterministic, no data.\n\n")

bs <- c(7.4, 14.9, 25, 40, 65, 100, 150, 250)
sd_h <- \(a) CL / (2 * sqrt(2 * a + 1))
tibble(b_shape = bs, stage_sd_h = sd_h(bs), cells = sd_h(bs) / (CL / 96),
       d_sd_per_10 = sd_h(bs) - sd_h(bs + 10)) |>
    mutate(across(-b_shape, \(x) round(x, 3))) |>
    print(n = Inf)

cat("\n=== 2. does the real posterior hit the upper bound? ===\n\n")
cat("Cells: pooled over the 14 grp_init groups and all post-warmup draws of\n")
cat("one real-data fit. `q50`/`q90`/`q975` are quantiles of b_shape,\n")
cat(sprintf("`p_gt_100` and `p_gt_200` the posterior probability above those\n"))
cat(sprintf("values, and `p_gt_240` the mass within 4%% of the max_shape = %d\n", MAX_SHAPE))
cat("bound -- if that is non-negligible the posterior is truncated and its\n")
cat("mean understates where the likelihood wants to go. Dimensionless.\n\n")

map(c("no_pool", "np_wide_bshape", "np_wide_both"), \(cfg) {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    m <- as.matrix(f, pars = "b_shape")
    rm(f); invisible(gc())
    v <- as.numeric(m)
    tibble(fit = cfg, mean = mean(v),
           q50 = quantile(v, 0.5), q90 = quantile(v, 0.9),
           q975 = quantile(v, 0.975),
           p_gt_100 = mean(v > 100), p_gt_200 = mean(v > 200),
           p_gt_240 = mean(v > 0.96 * MAX_SHAPE))
}) |> list_rbind() |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(width = Inf)

cat("\n  max_shape =", MAX_SHAPE, "\n")


cat("\n=== 3. what starting-stage age range does a given b_shape mean? ===\n\n")
cat("Inverts the mapping using the exact Beta quantiles rather than the\n")
cat("normal approximation. Cells: `age_range_h` is the width of the central\n")
cat("99% of the starting-stage distribution in hours -- the 0.5% to 99.5%\n")
cat("quantiles of Beta(b_shape, b_shape) scaled by cycle_length -- which is\n")
cat("the quantity mmcm.pdf's Figure 1 states as `a parasite age range\n")
cat("spanning 9 h (99% of the first wave of bursting)`. `b_shape` is the\n")
cat("value reproducing that range. Rows are cycle_length, to show how much\n")
cat("the answer depends on it. Deterministic, no data.\n\n")

b_for_range <- function(W, cl) {
    f <- \(a) (qbeta(0.995, a, a) - qbeta(0.005, a, a)) * cl - W
    if (f(2) < 0) return(NA_real_)
    uniroot(f, c(2, 1e4), tol = 1e-8)$root
}

expand_grid(age_range_h = c(6, 9, 12), cycle_length = c(43.7, 45.3, 46.1)) |>
    mutate(b_shape = map2_dbl(age_range_h, cycle_length, b_for_range)) |>
    mutate(b_shape = round(b_shape, 1)) |>
    pivot_wider(names_from = cycle_length, values_from = b_shape,
                names_prefix = "cl_") |>
    print(n = Inf)

cat("\nAnd the reverse, for the values in play:\n\n")
tibble(b_shape = c(14.9, 48.3, 65.4, 84, 100, 250)) |>
    mutate(age_range_h = (qbeta(0.995, b_shape, b_shape) -
                          qbeta(0.005, b_shape, b_shape)) * CL,
           age_range_h = round(age_range_h, 2)) |>
    print(n = Inf)
