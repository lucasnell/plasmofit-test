## Reader for the single-unit study (_scripts/unitsim.R): one row per merged
## fit, then per-chain cost. Descriptive for now: the study's reading rule is
## written into unitsim.R's header once the design is agreed, and this reader
## gains it then. First used 2026-10-21 on the timing probe (SLURM 32409).
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   UNITSIM_MODE=merge Rscript --vanilla _scripts/unitsim.R
##   Rscript --vanilla _scripts/unitsim-read.R > _data/unitsim-read-<label>-<date>.txt

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 140)
fs <- Sys.glob("_data/unitsim/sum-*.rds")
x <- bind_rows(lapply(fs, readRDS))
cat(length(fs), "merged fits\n\n")

cat("Cells: one row per fit (4 chains x 1000 post-warmup draws). max_rhat over\n",
    "all monitored quantities; true_cl, cl_mean (posterior mean), cl_sd\n",
    "(posterior sd), bias = cl_mean - true_cl, all hours, positive = long;\n",
    "in_80 = truth inside the central 80% posterior interval; bs_logratio =\n",
    "log(posterior mean b_shape / true b_shape); lp_apart = largest distance of\n",
    "a chain's mean lp__ from the median of the other three; div and\n",
    "treedepth-at-cap (12) in % of 4000 post-warmup iterations.\n\n", sep = "")
x |> transmute(name, unit, max_rhat = round(max_rhat, 3), worst, pass,
               true_cl = round(true_cl, 2), cl_mean = round(cl_mean, 2),
               cl_sd = round(cl_sd, 2), bias = round(bias, 3), in_80,
               bs_logratio = round(bs_logratio_gen, 3), lp_apart = round(lp_apart, 2),
               div_pct = round(div_pct, 2), td_cap_pct = round(treedepth_cap_pct, 2)) |>
    as.data.frame() |> print(row.names = FALSE)

cat("\nCells: per chain. hours = warmup + sampling; leapfrog steps per\n",
    "post-warmup iteration; ms per gradient = sampling seconds / sampling\n",
    "n_leapfrog; lp_mean and cl_mean = the chain's posterior means.\n\n", sep = "")
x |> select(name, per_chain) |> unnest(per_chain) |>
    mutate(across(c(lp_mean, cl_mean), \(v) round(v, 2)),
           across(c(hours, warmup_hours), \(v) round(v, 3)),
           leapfrog_per_iter = round(leapfrog_per_iter), ms_per_grad = round(ms_per_grad, 2)) |>
    as.data.frame() |> print(row.names = FALSE)
