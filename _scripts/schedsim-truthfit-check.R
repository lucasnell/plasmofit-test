# Option 2, per replicate: can this design tell b_shape = 15 from b_shape = 65?
#
# The by-arm tables in wockner-schedule-sim-analyze.R average over replicates.
# This prints the three SCHEDSIM_TRUTH_FIT=pl_wide_both replicates one at a
# time, against the same arm's replicates fitted to the pooled_cl truth, so
# the two truths can be read side by side under IDENTICAL priors (arm
# wide_nuis: sd_log_b_shape 1.5, sd_log10_total0 1).
#
# Reads the RES summaries plus each replicate's 55 MB fit, one at a time,
# because the RES object does not store per-group nuisance estimates. Minutes:
#     Rscript --vanilla _scripts/schedsim-truthfit-check.R \
#         2>&1 | tee _data/schedsim-truthfit-check.log

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({ library(tidyverse) })

## rep4-5 are included for the pooled_cl truth; only rep1-3 exist for the
## pl_wide_both truth, so the side-by-side comparison uses rep1-3 of each.
f <- list.files("_data", "^wock-schedsim-RES-wide_nuis-.*[.]rds$",
                full.names = TRUE)
res <- map(f, read_rds)

cat("Cells: one row per replicate of arm `wide_nuis`, which fits with\n")
cat("sd_log_b_shape = 1.5 and sd_log10_total0 = 1 in every row. `truth_fit`\n")
cat("is the saved fit the data were simulated from. `b_shape_true` and\n")
cat("`b_shape_est` are that parameter's value averaged over 14 grp_init\n")
cat("groups -- the truth used to simulate, and the posterior mean recovered;\n")
cat("`rel` is (est - true)/true, so 0 is perfect recovery. `cl_true` and\n")
cat("`cl_est` are the population cycle_length in hours and `cl_bias` is\n")
cat("est minus true. Not averaged over anything: one fit per row.\n\n")

map(res, \(r) {
    tibble(config = r$config,
           truth_fit = if (grepl("-truth", r$config)) sub(".*-truth", "", r$config)
                       else "pooled_cl",
           b_shape_true = mean(r$truth_nuisance$b_shape),
           b_shape_est = {
               f <- read_rds(sprintf("_data/wock-schedsim-fit-%s.rds", r$config))
               v <- mean(unname(colMeans(as.matrix(f, pars = "b_shape"))))
               rm(f); invisible(gc()); v
           },
           cl_true = r$true_cl,
           cl_est = mean(r$per_trial$cl_mean),
           max_rhat = r$max_rhat %||% NA_real_)
}) |>
    list_rbind() |>
    mutate(rel = (b_shape_est - b_shape_true) / b_shape_true,
           cl_bias = cl_est - cl_true) |>
    select(config, truth_fit, b_shape_true, b_shape_est, rel,
           cl_true, cl_est, cl_bias, max_rhat) |>
    arrange(truth_fit, config) |>
    mutate(across(where(is.numeric), \(x) round(x, 3))) |>
    print(n = Inf, width = Inf)
