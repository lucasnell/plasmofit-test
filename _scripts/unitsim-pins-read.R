## Reader for the single-unit fixed-pin grid (_scripts/unitsim.R header,
## "FIXED-PIN GRID"; design and reading rule fixed 2026-10-27 before any of
## its fits ran). Candidates per (window, rep) are the 11 fixed-b_shape fits
## of one dataset: pins truth / x2 / half at fit n_c 288, 384, 512, and the
## spread-matched pin at 288 and 512. The reference is pin truth at fit n_c
## 384, which is the main study's fix fit. elpd comes from the per-fit
## diag-<fit>.rds files of _scripts/unitsim-diag.R (compute mode).
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   UNITSIM_MODE=merge Rscript --vanilla _data/unitsim/unitsim-run-<n>.R
##   UNITSIM_DIAG=compute UNITSIM_DIAG_MATCH='-fix-rep[0-9]+(-pin|-f|$)' Rscript --vanilla _scripts/unitsim-diag.R
##   UNITSIM_PIN_PRIOR=<prod|truth|wide> Rscript --vanilla _scripts/unitsim-pins-read.R > _data/unitsim-read-pins-<date>.txt

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
options(width = 160)
DIR <- Sys.getenv("UNITSIM_DIR", "_data/unitsim/study/")
if (!endsWith(DIR, "/")) DIR <- paste0(DIR, "/")
PRIOR <- Sys.getenv("UNITSIM_PIN_PRIOR", "prod")
N_REP <- 60L; READABLE <- 45L
SD_TRUE <- 0.1276      # truth's total stage-spread sd at 168 h, cycles (unitsim-pinplan.R)

x <- bind_rows(lapply(Sys.glob(paste0(DIR, "sum-*.rds")), readRDS))
for (v in c("prior", "truth_set", "pin")) if (!v %in% names(x)) x[[v]] <- NA_character_
if (!"fit_nc" %in% names(x)) x$fit_nc <- NA_integer_
x <- x |>
    mutate(prior = coalesce(prior, "prod"), truth_set = coalesce(truth_set, "main"),
           pin = coalesce(pin, "truth"), fit_nc = coalesce(fit_nc, 384L)) |>
    filter(truth_set == "main", gen_nc == 384L, method == "fix", prior == PRIOR) |>
    mutate(cell = paste0(pin, "-f", fit_nc), k = 168 / true_cl)
## b_pin is missing from summaries merged before the grid (all truth pins).
if (!"b_pin" %in% names(x)) x$b_pin <- NA_real_
cat(nrow(x), "fixed-b_shape fits (gen n_c 384, main truth, prior", PRIOR, ") in", DIR, "\n\n")

## ---- 1. Per cell ----------------------------------------------------------
cat("1. PER CELL. n = merged fits; pass = max R-hat < 1.05; readable = >= ", READABLE,
    " of ", N_REP, " pass.\n",
    "bias = posterior mean cycle_length - truth (h), positive = long: mean and se over\n",
    "passing fits; cover80 = share of passing fits whose 80% interval holds the truth.\n",
    "spread_sd = the fit's total stage-spread sd at 168 h (cycles) from its pin and fit\n",
    "n_c; spread_diff = minus the truth's ", SD_TRUE, ". d_ref = paired bias minus the\n",
    "reference's (pin truth, fit 384), identical data, both passing: mean, se, pairs.\n\n", sep = "")
vi <- function(b) 1 / (4 * (2 * b + 1))
ref <- filter(x, pin == "truth", fit_nc == 384L, pass) |> select(window, rep, bias_ref = bias)
cells <- x |>
    left_join(ref, by = c("window", "rep")) |>
    group_by(window, fit_nc, pin) |>
    summarise(b_shape = first(na.omit(b_pin)), n = n(), n_pass = sum(pass),
              mean_bias = mean(bias[pass]), se = sd(bias[pass]) / sqrt(sum(pass)),
              cover80 = mean(in_80[pass]),
              d_ref = mean((bias - bias_ref)[pass & !is.na(bias_ref)]),
              d_ref_se = sd((bias - bias_ref)[pass & !is.na(bias_ref)]) /
                  sqrt(sum(pass & !is.na(bias_ref))),
              pairs = sum(pass & !is.na(bias_ref)), k = first(k), .groups = "drop") |>
    mutate(b_shape = if_else(is.na(b_shape) & pin == "truth", 20.098, b_shape),
           spread_sd = sqrt(vi(b_shape) + k / fit_nc), spread_diff = spread_sd - SD_TRUE,
           readable = n_pass >= READABLE) |>
    select(-k) |>
    mutate(across(where(is.double), \(v) round(v, 4)))
print(as.data.frame(cells), row.names = FALSE)

## ---- 2. Selection ----------------------------------------------------------
dfiles <- sprintf("%sdiag-%s.rds", DIR, x$name)
have <- file.exists(dfiles)
cat(sprintf("\n2. SELECTION. elpd available for %d of %d fits.\n", sum(have), nrow(x)))
if (sum(have) == 0L) quit(save = "no")
dg <- lapply(dfiles[have], readRDS); names(dg) <- x$name[have]
xs <- x[have, ] |> filter(pass) |>
    mutate(elpd = vapply(dg[name], `[[`, numeric(1), "elpd_loo"),
           ll_mean = vapply(dg[name], `[[`, numeric(1), "ll_mean"),
           absbias = abs(bias))
cat("Per (window, rep), candidates = the passing fits of that dataset (n_cand).\n",
    "selected = highest elpd_loo; oracle = smallest |bias| among candidates;\n",
    "n_close = candidates within 2 se (pointwise difference) of the selected,\n",
    "including it (1 = decisive); rho = Spearman correlation over candidates of\n",
    "elpd with |bias| (negative = higher elpd goes with less bias).\n\n", sep = "")
sel <- xs |> group_by(window, rep) |> group_modify(\(g, k) {
    i <- which.max(g$elpd); best <- g$name[i]
    se_d <- vapply(g$name, \(nm) if (nm == best) 0 else
        sd(dg[[nm]]$pointwise - dg[[best]]$pointwise) * sqrt(length(dg[[best]]$pointwise)),
        numeric(1))
    tibble(n_cand = nrow(g), selected = g$cell[i], sel_absbias = g$absbias[i],
           ref_in = any(g$cell == "truth-f384"),
           ref_absbias = if (any(g$cell == "truth-f384")) g$absbias[g$cell == "truth-f384"] else NA_real_,
           oracle = g$cell[which.min(g$absbias)], oracle_absbias = min(g$absbias),
           n_close = sum(g$elpd >= g$elpd[i] - 2 * se_d),
           rho = suppressWarnings(cor(g$elpd, g$absbias, method = "spearman")),
           sel_ll = g$cell[which.max(g$ll_mean)],
           sel_ll_absbias = g$absbias[which.max(g$ll_mean)],
           rho_ll = suppressWarnings(cor(g$ll_mean, g$absbias, method = "spearman")))
}) |> ungroup()

cat("Summary per window (means over reps):\n")
sel |> group_by(window) |>
    summarise(reps = n(), n_cand = mean(n_cand), ref_selected = sum(selected == "truth-f384"),
              oracle_selected = sum(selected == oracle),
              sel_absbias = mean(sel_absbias), ref_absbias = mean(ref_absbias, na.rm = TRUE),
              oracle_absbias = mean(oracle_absbias), n_close = mean(n_close),
              rho = mean(rho, na.rm = TRUE), .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)

cat("\nHow often each cell is selected (by elpd) and is the oracle, per window:\n")
full_join(count(sel, window, cell = selected, name = "selected"),
          count(sel, window, cell = oracle, name = "oracle"), by = c("window", "cell")) |>
    mutate(across(c(selected, oracle), \(v) coalesce(v, 0L))) |> arrange(window, desc(selected)) |>
    as.data.frame() |> print(row.names = FALSE)

cat("\nDisplay only: the same selection by posterior mean summed log-likelihood (in-sample).\n")
sel |> group_by(window) |>
    summarise(reps = n(), ref_selected = sum(sel_ll == "truth-f384"),
              sel_absbias = mean(sel_ll_absbias), rho = mean(rho_ll, na.rm = TRUE),
              .groups = "drop") |>
    mutate(across(where(is.double), \(v) round(v, 3))) |> as.data.frame() |> print(row.names = FALSE)
count(sel, window, cell = sel_ll, name = "selected") |> arrange(window, desc(selected)) |>
    as.data.frame() |> print(row.names = FALSE)
