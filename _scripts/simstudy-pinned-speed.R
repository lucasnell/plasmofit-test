## Speed of the simulation study's pinned chains (jobs 31864/31865, task IDs
## in _data/simstudy/pinned-tasks.txt) against the packed fit n_c 384 chains
## of SLURM 30875, from the chain files (descriptive; added 2026-10-19).
## Gradients = n_leapfrog summed over warmup and sampling; time =
## rstan::get_elapsed_time (warmup + sample). Also projects the end of the
## pinned queue: remaining chains spread over 12 workers at the median
## pinned chain time.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/simstudy-pinned-speed.R > _data/simstudy-pinned-speed-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages({library(rstan); library(dplyr)})
options(width = 120)
pinned <- as.integer(readLines("_data/simstudy/pinned-tasks.txt"))
FITS_PER_REP <- 48L; N_CHAINS <- 8L
cls <- function(id) {
    f <- (id - 1L) %/% N_CHAINS + 1L; j <- (f - 1L) %% FITS_PER_REP + 1L
    ifelse(j <= 24L, 384L, 192L)
}
## task ID -> chain file, via the log line the frozen script writes
logf <- sprintf("_data/simstudy/logs/%d.out", 1:1920)
rows <- list()
for (id in which(file.exists(logf))) {
    if (cls(id) != 384L) next
    L <- readLines(logf[id], warn = FALSE)
    w <- grep("^wrote ", L, value = TRUE)
    if (length(w) == 0) next
    out <- trimws(sub("^wrote ", "", w[1]))
    if (!file.exists(out)) next   # merged runs may have removed it
    ch <- readRDS(out)
    sp <- get_sampler_params(ch$fit, inc_warmup = TRUE)[[1]]
    rows[[length(rows) + 1]] <- data.frame(
        id = id, layout = ifelse(id %in% pinned, "pinned", "packed"),
        method = ifelse(grepl("-wide-", out), "wide", "fixed"),
        grads = sum(sp[, "n_leapfrog__"]), secs = sum(get_elapsed_time(ch$fit)))
}
d <- bind_rows(rows)

cat("Cells: fit n_c 384 chains with a finished chain file, per chain, median\n")
cat("(range) over chains; gradients in thousands; ms per gradient = chain\n")
cat("seconds / chain gradients. Unpaired: different simulated datasets.\n\n")
tab <- d |> group_by(layout, method) |>
    summarise(chains = n(),
              k_grads = sprintf("%.0f (%.0f-%.0f)", median(grads) / 1e3, min(grads) / 1e3, max(grads) / 1e3),
              ms_per_grad = sprintf("%.1f (%.1f-%.1f)", median(1e3 * secs / grads),
                                    min(1e3 * secs / grads), max(1e3 * secs / grads)),
              hours = sprintf("%.2f (%.2f-%.2f)", median(secs) / 3600, min(secs) / 3600, max(secs) / 3600),
              .groups = "drop")
print(as.data.frame(tab), row.names = FALSE)

p <- d |> filter(layout == "pinned")
done <- nrow(p); left <- length(pinned) - done
h <- median(p$secs) / 3600
cat(sprintf("\npinned: %d of %d chains finished, %d left (running or unclaimed)\n",
            done, length(pinned), left))
cat(sprintf("crude projection: %d left / 12 workers x %.2f h median = %.1f h of pinned work from now\n",
            left, h, left / 12 * h))
cat(sprintf("now %s\n", format(Sys.time(), "%Y-%m-%d %H:%M")))
