## Where the simulation study's gate failures sit (descriptive; added
## 2026-10-18 after the first partial read, so NOT part of the pre-registered
## reading rule in simstudy.R). Bins of max R-hat are descriptive labels,
## chosen after seeing the data: "marginal" 1.05-1.15, "severe" >= 1.15.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/simstudy-health.R > _data/simstudy-health-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(dplyr); library(readr)})
options(width = 120)
S <- bind_rows(lapply(Sys.glob("_data/simstudy/sum-*.rds"), read_rds))
cat("fits merged:", nrow(S), "\n\n")

cat("Cells: counts of merged fits by max R-hat over all quantities (8 chains x 500\n")
cat("draws); lp_only = failing fits whose ONLY quantity at or above 1.05 is lp__;\n")
cat("med_lp_spread = median over fits of the range of per-chain mean lp__.\n\n")
tab <- S |>
    mutate(bin = case_when(pass ~ "pass", max_rhat < 1.15 ~ "marginal", TRUE ~ "severe"),
           lp_only = !pass & worst == "lp__" & n_rhat_over == 1) |>
    group_by(method, fit_nc) |>
    summarise(n = n(), pass = sum(bin == "pass"), marginal = sum(bin == "marginal"),
              severe = sum(bin == "severe"), lp_only = sum(lp_only),
              med_lp_spread = round(median(lp_spread), 1),
              med_div_pct = round(median(div_pct), 2),
              med_hours = round(median(hours_max_chain), 2), .groups = "drop")
print(as.data.frame(tab), row.names = FALSE)

cat("\nSame bins by cell (truth, gen_nc, fit_nc, bs_scale, method), all replicates merged:\n\n")
cell <- S |>
    group_by(truth, gen_nc, fit_nc, bs_scale, method) |>
    summarise(n = n(), pass = sum(pass),
              severe = sum(!pass & max_rhat >= 1.15),
              max_rhats = paste(sprintf("%.2f", sort(max_rhat)), collapse = " "),
              .groups = "drop")
print(as.data.frame(cell), row.names = FALSE)

cat("\nTrial most often worst among severe failures (cycle_length/eta_cl/logit_cl index):\n")
sev <- S |> filter(!pass, max_rhat >= 1.15)
idx <- suppressWarnings(as.integer(sub(".*\\[(\\d+).*", "\\1",
                                       sev$worst[grepl("cycle_length|eta_cl|logit_cl", sev$worst)])))
trials <- S$per_trial[[1]]$trial
print(sort(table(trials[idx]), decreasing = TRUE))
