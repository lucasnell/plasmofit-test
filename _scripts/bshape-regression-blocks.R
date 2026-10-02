# Where does the b_shape regression's excess sit?
#
# The test (old build vs new build, same seed) put 17.5% of entries beyond 2
# combined MCSE. The null (new build vs new build, seed alone) puts 4.1%
# there. Both have b_shape free and heavy-tailed, so "the posterior mean of a
# heavy-tailed parameter is a high-variance statistic" cannot by itself
# explain the gap -- the null controls for exactly that and does not
# reproduce it.
#
# What remains between them is the BUILD. This splits the comparison by
# parameter block to ask whether the excess is concentrated where the code
# change touches (b_shape, and whatever correlates with it) or spread evenly,
# which would point at generic rebuild drift instead.
#
# Reads three 55 MB fits one at a time.
#     srun -N 1 -n 1 -c 4 --mem=48G Rscript --vanilla \
#         _scripts/bshape-regression-blocks.R 2>&1 | tee _data/bshape-regression-blocks.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

get1 <- function(cfg) {
    f <- read_rds(sprintf("_data/wock-fit-%s.rds", cfg))
    ss <- rstan::summary(f)$summary
    rm(f); invisible(gc())
    tibble(par = rownames(ss), mean = ss[, "mean"], se = ss[, "se_mean"])
}

pair <- function(a, b, label) {
    inner_join(get1(a), get1(b), by = "par", suffix = c("_n", "_o")) |>
        filter(par != "lp__", !startsWith(par, "log_lik"),
               is.finite(se_n), is.finite(se_o), se_n > 0, se_o > 0) |>
        mutate(z = (mean_n - mean_o) / sqrt(se_n^2 + se_o^2),
               block = sub("\\[.*$", "", par),
               cmp = label)
}

cat("=== reading ===\n")
## testB is the same build contrast with b_shape PINNED by a tight prior
## instead of free, which is the one structural difference between the two
## test pairs.
d <- bind_rows(
    pair("np_wide_total0-rebuild", "np_wide_total0",         "test"),
    pair("np_wide_total0-null",    "np_wide_total0-rebuild", "null"),
    pair("np_bs400-rebuild",       "np_bs400",               "testB"))

cat("\nCells: one row per parameter block, over the scalar entries both fits\n")
cat("of a pair share. `z` is the difference in posterior means over the two\n")
cat("runs' combined Monte Carlo standard error. `test` is old build vs new\n")
cat("build at the SAME seed; `null` is new build vs new build at a different\n")
cat("seed, so it is run-to-run variation alone. `mean|z|` averages |z| in the\n")
cat("block and `f>2` is the share beyond 2. `ratio` is the test's mean|z|\n")
cat("over the null's: near 1 means the build adds nothing beyond the seed.\n")
cat("Dimensionless. n is per block, identical in both comparisons.\n\n")

d |>
    summarise(.by = c(cmp, block), n = n(), mabs = mean(abs(z)),
              mx = max(abs(z)), f2 = mean(abs(z) > 2)) |>
    pivot_wider(names_from = cmp, values_from = c(n, mabs, mx, f2)) |>
    transmute(block, n = n_null,
              mabs_null, mabs_test, mabs_testB,
              ratio_A = mabs_test / mabs_null,
              ratio_B = mabs_testB / mabs_null) |>
    arrange(desc(ratio_A)) |>
    mutate(across(where(is.numeric), \(x) round(x, 2))) |>
    print(n = Inf, width = Inf)

cat("\nOverall:\n")
d |> summarise(.by = cmp, n = n(), mean_abs_z = mean(abs(z)),
               max_abs_z = max(abs(z)), frac_gt2 = mean(abs(z) > 2)) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |> print(width = Inf)
