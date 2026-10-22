## Reader for the single-unit study (_scripts/unitsim.R; design and reading
## rule in its header, fixed 2026-10-22 before any study fit ran). Reads the
## merged summaries in UNITSIM_DIR (default _data/unitsim/study/). Cells are
## keyed by `key` = window (w168 / w192) for the study; the 2026-10-21 probes
## (UNITSIM_DIR=_data/unitsim) predate windows and are keyed by unit, and for
## them the per-fit and per-chain tables are the useful part.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   UNITSIM_MODE=merge Rscript --vanilla _scripts/unitsim.R
##   Rscript --vanilla _scripts/unitsim-read.R > _data/unitsim-read-<label>-<date>.txt

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 150)
DIR <- Sys.getenv("UNITSIM_DIR", "_data/unitsim/study/")
N_REP <- 60L; READABLE <- 45L
fs <- Sys.glob(file.path(DIR, "sum-*.rds"))
x <- bind_rows(lapply(fs, readRDS))
if (!"window" %in% names(x)) x$window <- NA_integer_
x <- x |>
    mutate(key = if_else(is.na(window), unit, paste0("w", window)),
           rhat_bin = case_when(pass ~ "pass", max_rhat < 1.15 ~ "marginal", TRUE ~ "severe"),
           apart = lp_apart > 3)
cat(length(fs), "merged fits in", DIR, "\n\n")

## ---- 1. Per cell ----------------------------------------------------------
cat("1. PER CELL. Cells: n = merged fits; pass = max R-hat < 1.05 (4 chains);\n",
    "marginal 1.05-1.15, severe >= 1.15; apart = failing fits with a chain's\n",
    "mean lp__ > 3 from the median of the other three. bias = posterior mean\n",
    "cycle_length - truth, hours, positive = long: mean and se (sd / sqrt(n))\n",
    "over PASSING fits, unpaired across cells. readable = >= ", READABLE, " of ", N_REP,
    " pass.\n", "cover80 = share of passing fits whose central 80% interval holds\n",
    "the truth. bs_log = est arm, mean log(posterior mean b_shape / truth).\n",
    "all_bias = mean over ALL fits, display only.\n\n", sep = "")
cells <- x |> group_by(key, gen_nc, method) |>
    summarise(n = n(), pass = sum(pass), marginal = sum(rhat_bin == "marginal"),
              severe = sum(rhat_bin == "severe"), apart = sum(!pass & apart),
              bias = mean(bias[pass]), se = sd(bias[pass]) / sqrt(sum(pass)),
              cover80 = mean(in_80[pass]),
              bs_log = if (first(method) == "est") mean(bs_logratio_gen[pass]) else NA_real_,
              all_bias = mean(bias), .groups = "drop") |>
    mutate(readable = pass >= READABLE,
           across(c(bias, se, cover80, bs_log, all_bias), \(v) round(v, 3)))
print(as.data.frame(cells), row.names = FALSE)

cat("\nPassing against failing fits, display only. Cells: mean bias (hours) and\n",
    "count, by cell.\n\n", sep = "")
x |> group_by(key, gen_nc, method, pass) |>
    summarise(n = n(), mean_bias = round(mean(bias), 3), .groups = "drop") |>
    pivot_wider(names_from = pass, values_from = c(n, mean_bias), names_prefix = "pass_") |>
    as.data.frame() |> print(row.names = FALSE)

## ---- 2. Paired contrasts ---------------------------------------------------
cat("\n2. PAIRED CONTRASTS. Cells: per replicate, difference of posterior mean\n",
    "cycle_length (hours) between two fits of IDENTICAL noise draws, both\n",
    "passing; mean, se, range over pairs, and how many are negative.\n",
    "w192 - w168: same truth and the same noise on the 34 shared\n",
    "observations; the effect of extending the window to 192 h.\n",
    "est - fix: same data, b_shape estimated minus fixed at truth (cost of\n",
    "estimating it). gen 1024 - 384: same noise, different noiseless\n",
    "trajectories (partially paired).\n\n", sep = "")
pair_sum <- function(d) d |>
    summarise(pairs = n(), mean = round(mean(diff), 3), se = round(sd(diff) / sqrt(n()), 3),
              lo = round(min(diff), 3), hi = round(max(diff), 3), negative = sum(diff < 0),
              .groups = "drop")
xp <- filter(x, pass)
w_diff <- inner_join(filter(xp, key == "w192"), filter(xp, key == "w168"),
                     by = c("gen_nc", "method", "rep"), suffix = c("_l", "_s")) |>
    mutate(diff = cl_mean_l - cl_mean_s) |> group_by(gen_nc, method) |> pair_sum()
cat("w192 - w168\n"); print(as.data.frame(w_diff), row.names = FALSE)
est_fix <- inner_join(filter(xp, method == "est"), filter(xp, method == "fix"),
                      by = c("key", "gen_nc", "rep"), suffix = c("_e", "_f")) |>
    mutate(diff = cl_mean_e - cl_mean_f) |> group_by(key, gen_nc) |> pair_sum()
cat("\nest - fix\n"); print(as.data.frame(est_fix), row.names = FALSE)
g_diff <- inner_join(filter(xp, gen_nc == 1024L), filter(xp, gen_nc == 384L),
                     by = c("key", "method", "rep"), suffix = c("_h", "_l")) |>
    mutate(diff = bias_h - bias_l) |> group_by(key, method) |> pair_sum()
cat("\ngen 1024 - 384 (difference of bias, since the truths match)\n")
print(as.data.frame(g_diff), row.names = FALSE)

## ---- 3. Per fit and per chain (probes; or to inspect) ------------------------
if (nrow(x) <= 24L || nzchar(Sys.getenv("UNITSIM_DETAIL"))) {
    cat("\n3. PER FIT. max_rhat over all monitored quantities; cl_mean, cl_sd =\n",
        "posterior mean and sd of cycle_length; bias = cl_mean - truth (h);\n",
        "lp_apart as above; div and treedepth-at-cap (12) in % of post-warmup\n",
        "iterations.\n\n", sep = "")
    x |> transmute(name, unit, max_rhat = round(max_rhat, 3), worst, pass,
                   true_cl = round(true_cl, 2), cl_mean = round(cl_mean, 2),
                   cl_sd = round(cl_sd, 2), bias = round(bias, 3), in_80,
                   bs_logratio = round(bs_logratio_gen, 3), lp_apart = round(lp_apart, 2),
                   div_pct = round(div_pct, 2), td_cap_pct = round(treedepth_cap_pct, 2)) |>
        as.data.frame() |> print(row.names = FALSE)
    cat("\nPer chain: hours = warmup + sampling; leapfrog steps per post-warmup\n",
        "iteration; ms per gradient = sampling seconds / sampling n_leapfrog.\n\n", sep = "")
    x |> select(name, per_chain) |> unnest(per_chain) |>
        mutate(across(c(lp_mean, cl_mean), \(v) round(v, 2)),
               across(c(hours, warmup_hours), \(v) round(v, 3)),
               leapfrog_per_iter = round(leapfrog_per_iter), ms_per_grad = round(ms_per_grad, 2)) |>
        as.data.frame() |> print(row.names = FALSE)
}
