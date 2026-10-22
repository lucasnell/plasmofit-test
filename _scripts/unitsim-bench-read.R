## Reads the layout benchmark for single-unit chains (2026-10-22; jobs in
## _data/unitsim/bench-jobs.txt): the same 64 chain tasks (the single-unit
## study's first 64 tokens: reps 1-2 of all 8 cells; fix arm 2000 + 2000
## iterations, est arm 700 + 1000) run twice on 128 CPUs = 16 L3 groups:
##   w4  one chain per L3 group, 4 physical cores each (16 workers; the
##       always-pin layout);
##   w2  two chains per L3 group, 2 physical cores each (32 workers).
## Lucas decides the layout from these numbers; no rule is applied here.
## The runs are sequential, so other users' load on the node may differ.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/unitsim-bench-read.R > _data/unitsim-bench-read-<date>.txt

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(rstan)})
options(width = 140)
N_GROUPS <- 16L

chains <- function(layout) {
    fs <- Sys.glob(sprintf("_data/unitsim/bench-%s/chain-*-c*.rds", layout))
    bind_rows(lapply(fs, \(f) {
        x <- readRDS(f); el <- get_elapsed_time(x$fit)
        sp <- get_sampler_params(x$fit, inc_warmup = FALSE)[[1]]   # no warmup saved
        tibble(layout = layout,
               token = sub("^chain-(.*)\\.rds$", "\\1", basename(f)),
               method = x$cf$method, unit = x$truth$unit,
               hours = sum(el) / 3600,
               ms_per_grad = 1e3 * el[, "sample"] / sum(sp[, "n_leapfrog__"]),
               leapfrog_per_iter = mean(sp[, "n_leapfrog__"]))
    }))
}
x <- bind_rows(chains("w4"), chains("w2"))

## Wall time of each run: first task start to last worker finish, from the
## driver log.
span <- function(layout) {
    L <- readLines(Sys.glob(sprintf("_data/unitsim/bench-%s/logs/pin-*.out", layout))[1])
    ts <- as.POSIXct(substr(grep("^20[0-9-]+ [0-9:]+ ", L, value = TRUE), 1, 19))
    as.numeric(difftime(max(ts), min(ts), units = "hours"))
}
wall <- c(w4 = span("w4"), w2 = span("w2"))

cat("Chains finished:", table(x$layout)[c("w4", "w2")], "(of 64 each)\n\n")
cat("Cells: per layout, all 64 chains; wall = first start to last finish (h);\n",
    "throughput = chains per L3 group per hour = 64 / (16 groups x wall);\n",
    "chain hours and ms per gradient (sampling seconds / sampling n_leapfrog):\n",
    "median (range) over chains.\n\n", sep = "")
mr <- function(v, d = 3) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), median(v), min(v), max(v))
x |> group_by(layout) |>
    summarise(chains = n(), chain_hours = mr(hours), ms_per_grad = mr(ms_per_grad, 2),
              leapfrog_per_iter = mr(leapfrog_per_iter, 0), .groups = "drop") |>
    mutate(wall_h = round(wall[layout], 3),
           throughput = round(64 / (N_GROUPS * wall[layout]), 2)) |>
    arrange(desc(layout)) |> as.data.frame() |> print(row.names = FALSE)

cat("\nPaired by chain (same data, seed, and settings in both layouts). Cells:\n",
    "ratio w2 / w4 of chain hours and of ms per gradient, median (range) over\n",
    "chains, by arm. A chain-hours ratio r means w2 has 2 / r times w4's\n",
    "throughput per L3 group, ignoring tail effects.\n\n", sep = "")
x |> select(layout, token, method, unit, hours, ms_per_grad, leapfrog_per_iter) |>
    pivot_wider(names_from = layout, values_from = c(hours, ms_per_grad, leapfrog_per_iter)) |>
    filter(!is.na(hours_w4), !is.na(hours_w2)) |>
    group_by(method) |>
    summarise(pairs = n(), hours_ratio = mr(hours_w2 / hours_w4),
              ms_per_grad_ratio = mr(ms_per_grad_w2 / ms_per_grad_w4),
              leapfrog_ratio = mr(leapfrog_per_iter_w2 / leapfrog_per_iter_w4),
              .groups = "drop") |>
    as.data.frame() |> print(row.names = FALSE)
