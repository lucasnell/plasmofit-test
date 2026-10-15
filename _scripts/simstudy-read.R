## Read the cycle_length simulation study (_scripts/simstudy.R). The reading
## rule is fixed in that script's header, before any fit ran; this applies it.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/simstudy-read.R > _data/simstudy-read-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})
options(width = 120)

## SIMSTUDY_SUMGLOB lets a smoke test read its own summaries.
fs <- Sys.glob(Sys.getenv("SIMSTUDY_SUMGLOB", "_data/simstudy/sum-*.rds"))
cat("fits merged:", length(fs), "of 240\n\n")
if (length(fs) == 0) quit(save = "no")
S <- bind_rows(lapply(fs, read_rds))

cat("=== sampler health over all merged fits ===\n")
cat(sprintf("pass the gate (max R-hat < 1.05, 8 chains): %d of %d\n", sum(S$pass), nrow(S)))
cat(sprintf("divergent %% of post-warmup transitions: median %.2f, max %.2f\n",
            median(S$div_pct), max(S$div_pct)))
cat(sprintf("longest chain, hours: median %.2f, max %.2f; chain-hours in total %.0f\n\n",
            median(S$hours_max_chain), max(S$hours_max_chain), sum(S$hours_sum_chains)))

cat("=== bias by cell ===\n")
cat("Cells: bias = mean over 13 trials of (posterior mean cycle_length - truth),\n")
cat("hours, positive = estimate runs long; mean and range over the replicates\n")
cat("that pass the gate, n_pass of the replicates merged. NOT READABLE if fewer\n")
cat("than 3 replicates pass. Truth A: 45.012 h in every trial; truth B: per-trial\n")
cat("truths, mean 41.82 h. Unpaired across cells.\n\n")
cells <- S |>
    group_by(truth, gen_nc, fit_nc, bs_scale, method) |>
    summarise(n = n(), n_pass = sum(pass),
              bias = if (sum(pass) >= 3) mean(bias[pass]) else NA_real_,
              lo = if (sum(pass) >= 3) min(bias[pass]) else NA_real_,
              hi = if (sum(pass) >= 3) max(bias[pass]) else NA_real_,
              bs_logratio = if (sum(pass) >= 3) mean(bs_logratio[pass]) else NA_real_,
              .groups = "drop") |>
    mutate(across(c(bias, lo, hi, bs_logratio), \(x) round(x, 3)),
           note = ifelse(n_pass < 3, "NOT READABLE", ""))
print(as.data.frame(cells), row.names = FALSE)
cat("\nbs_logratio: wide arm only, mean over 14 units of log(posterior mean b_shape /\n")
cat("true b_shape); 0 = recovered, negative = estimated too low (less synchronous).\n\n")

pair_on <- c("truth", "gen_nc", "bs_scale", "rep")
## A partial read may hold one side of a pair only; add the missing column.
fill_cols <- function(x, cols) { for (cn in setdiff(cols, names(x))) x[[cn]] <- NA_real_; x }
cat("=== paired: wide - fixed, same simulated data and fit_nc ===\n")
cat("Cells: mean (range) over replicates where both fits pass, hours.\n\n")
pw <- S |> filter(pass) |>
    select(all_of(pair_on), fit_nc, method, bias) |>
    pivot_wider(names_from = method, values_from = bias) |>
    fill_cols(c("wide", "fixed")) |>
    filter(!is.na(wide), !is.na(fixed)) |>
    mutate(d = wide - fixed) |>
    group_by(truth, gen_nc, fit_nc, bs_scale) |>
    summarise(n_pair = n(), mean = round(mean(d), 3),
              lo = round(min(d), 3), hi = round(max(d), 3), .groups = "drop")
print(as.data.frame(pw), row.names = FALSE)

cat("\n=== paired: fit_nc 384 - fit_nc 192, same simulated data and method ===\n")
cat("Cells: mean (range) over replicates where both fits pass, hours. At gen_nc 192\n")
cat("this is over-specifying n_c; at gen_nc 384, the correct fit minus the\n")
cat("under-specified one.\n\n")
pn <- S |> filter(pass) |>
    select(all_of(pair_on), method, fit_nc, bias) |>
    pivot_wider(names_from = fit_nc, values_from = bias, names_prefix = "f") |>
    fill_cols(c("f384", "f192")) |>
    filter(!is.na(f384), !is.na(f192)) |>
    mutate(d = f384 - f192) |>
    group_by(truth, gen_nc, bs_scale, method) |>
    summarise(n_pair = n(), mean = round(mean(d), 3),
              lo = round(min(d), 3), hi = round(max(d), 3), .groups = "drop")
print(as.data.frame(pn), row.names = FALSE)

cat("\n=== fits that fail the gate ===\n")
print(as.data.frame(S |> filter(!pass) |>
                        select(name, max_rhat, worst, n_rhat_over, lp_spread) |>
                        mutate(across(where(is.numeric), \(x) round(x, 3))) |>
                        arrange(name)), row.names = FALSE)
