## Split the simulation study's fit_nc 384 vs 192 run-time ratio into cost per
## gradient and number of gradients, from the merged fits (descriptive; added
## 2026-10-18 to answer whether a cheaper forward map at the n_c 384
## dispersion would shorten fits). Gradients = n_leapfrog summed over warmup
## and sampling; time = rstan::get_elapsed_time (warmup + sample), per chain.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/simstudy-cost.R > _data/simstudy-cost-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(dplyr)})
options(width = 120)
files <- Sys.glob("_data/simstudy/fit-*.rds")
rows <- lapply(files, function(f) {
    fit <- readRDS(f)
    sp <- get_sampler_params(fit, inc_warmup = TRUE)
    el <- get_elapsed_time(fit)
    nm <- sub("^fit-(.*)\\.rds$", "\\1", basename(f))
    p <- strsplit(nm, "-")[[1]]
    data.frame(name = nm, fit_nc = as.integer(sub("f", "", p[3])), method = p[5],
               chain = seq_along(sp),
               grads = vapply(sp, function(s) sum(s[, "n_leapfrog__"]), numeric(1)),
               secs = rowSums(el))
})
d <- bind_rows(rows)
pass <- bind_rows(lapply(Sys.glob("_data/simstudy/sum-*.rds"), readRDS)) |>
    transmute(name = sub("^sum-", "", name), pass)
d <- left_join(d, pass, by = "name")

cat("fits read:", length(files), "; chains:", nrow(d), "\n\n")
cat("Cells: per chain, median over chains in the group; gradients in thousands\n")
cat("(warmup 700 + sampling 500 iterations); ms per gradient = chain seconds /\n")
cat("chain gradients. Unpaired; all fits merged so far, passing or not.\n\n")
tab <- d |>
    group_by(fit_nc, method) |>
    summarise(chains = n(), k_grads = round(median(grads) / 1000, 1),
              ms_per_grad = round(median(1000 * secs / grads), 2),
              hours = round(median(secs) / 3600, 2), .groups = "drop")
print(as.data.frame(tab), row.names = FALSE)

cat("\nRatio fit_nc 384 / 192 of the medians above, by method:\n")
r <- tab |> group_by(method) |> arrange(fit_nc, .by_group = TRUE) |>
    summarise(hours = round(last(hours) / first(hours), 2),
              grads = round(last(k_grads) / first(k_grads), 2),
              ms_per_grad = round(last(ms_per_grad) / first(ms_per_grad), 2))
print(as.data.frame(r), row.names = FALSE)

cat("\nSame, restricted to fits that pass the gate:\n")
tp <- d |> filter(pass) |> group_by(fit_nc, method) |>
    summarise(chains = n(), k_grads = round(median(grads) / 1000, 1),
              ms_per_grad = round(median(1000 * secs / grads), 2),
              hours = round(median(secs) / 3600, 2), .groups = "drop")
print(as.data.frame(tp), row.names = FALSE)
