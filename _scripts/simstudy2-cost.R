## Cost per chain of the configuration the second simulation study fits
## (`_scripts/simstudy2.R`): no_pool, n_c 384, b_shape estimated under
## lognormal(2, 0.5), sd_log10_total0 1, adapt_delta 0.95, max_treedepth 12,
## 700 warmup + 1000 draws. No simulated fit has run under it, so two proxies
## (descriptive; added 2026-10-21 to size the study's pilot):
##   1. the real-data fits with exactly this configuration, pinned layout
##      (seed probe, SLURM 30858: np_wide_total0_nc384_t, seeds 20261013-15),
##      and its original unpinned run (SLURM 30842);
##   2. the old simulation study's fits of the same truth at fit n_c 384 with
##      b_shape prior sd 1.5 and the default sampler (truth B, b_shape x 1,
##      wide arm; 700 warmup + 500 draws).
## For contrast, the pinned real-data fits with b_shape pinned at 15 in every
## unit (default sampler). Gradients = n_leapfrog; time = get_elapsed_time.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/simstudy2-cost.R > _data/simstudy2-cost-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(dplyr)})
options(width = 140)

chain_row <- function(fit, group, file) {
    el <- get_elapsed_time(fit)
    ## Warmup draws are not saved in these files, so sampler parameters cover
    ## the post-warmup iterations only; cost per gradient is measured there.
    smp <- get_sampler_params(fit, inc_warmup = FALSE)[[1]]
    tibble(group = group, file = basename(file),
           hours = sum(el) / 3600,
           warmup_hours = el[, "warmup"] / 3600,
           ms_per_grad = 1e3 * el[, "sample"] / sum(smp[, "n_leapfrog__"]),
           leapfrog_per_iter = mean(smp[, "n_leapfrog__"]),
           treedepth_med = median(smp[, "treedepth__"]),
           treedepth_max = max(smp[, "treedepth__"]))
}

## wockner-fit.R chain files hold a list with $fit; older ones a bare stanfit.
as_fit <- function(x) if (inherits(x, "stanfit")) x else x$fit

rows <- list()
for (f in Sys.glob("_data/wock-chain-np_wide_total0_nc384_t-pinned-s*.rds"))
    rows[[f]] <- chain_row(as_fit(readRDS(f)), "real, b_shape estimated sd 0.5, tuned, pinned", f)
for (f in Sys.glob("_data/wock-chain-np_bs15_nc384-pinned*.rds"))
    rows[[f]] <- chain_row(as_fit(readRDS(f)), "real, b_shape pinned 15, default, pinned", f)
for (f in Sys.glob("_data/simstudy/chain-B-g*-f384-bs1.0-wide-rep*-c*.rds"))
    rows[[f]] <- chain_row(readRDS(f)$fit, "sim, b_shape estimated sd 1.5, default", f)
x <- bind_rows(rows)

## Original unpinned run: one merged stanfit, split by chain.
f0 <- readRDS("_data/wock-fit-np_wide_total0_nc384_t.rds")
el0 <- get_elapsed_time(f0)
sp0 <- get_sampler_params(f0, inc_warmup = FALSE)
orig <- tibble(chain = seq_along(sp0), hours = rowSums(el0) / 3600,
               leapfrog_per_iter = sapply(sp0, \(s) mean(s[, "n_leapfrog__"])),
               treedepth_med = sapply(sp0, \(s) median(s[, "treedepth__"])),
               treedepth_max = sapply(sp0, \(s) max(s[, "treedepth__"])))

cat("Cells: per chain, median (range) over chains. hours = warmup + sampling;\n",
    "warmup_hours = warmup alone; ms per gradient = sampling seconds / sampling\n",
    "n_leapfrog (warmup draws are not saved); leapfrog per iteration and\n",
    "treedepth over post-warmup iterations. Unpaired. Real-data chains: 700 +\n",
    "1000 iterations; simulated: 700 + 500.\n\n", sep = "")
mr <- function(v, d = 1) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"),
                                median(v), min(v), max(v))
x |> group_by(group) |>
    summarise(chains = n(), hours = mr(hours, 2), warmup_hours = mr(warmup_hours, 2),
              ms_per_grad = mr(ms_per_grad), leapfrog_per_iter = mr(leapfrog_per_iter, 0),
              treedepth = sprintf("med %g, max %g", median(treedepth_med), max(treedepth_max))) |>
    as.data.frame() |> print(row.names = FALSE)

cat("\nOriginal unpinned run of the real-data fit (SLURM 30842), per chain:\n\n")
print(as.data.frame(orig |> mutate(across(c(hours, leapfrog_per_iter), \(v) round(v, 1)))),
      row.names = FALSE)

cat("\nPer chain, pinned real-data fits with b_shape estimated:\n\n")
x |> filter(startsWith(group, "real, b_shape estimated")) |>
    select(-group) |> mutate(across(where(is.double), \(v) round(v, 1))) |>
    as.data.frame() |> print(row.names = FALSE)
