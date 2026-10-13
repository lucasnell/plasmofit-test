## Where does the time go in the n_c = 384 fits? Cost per chain =
## (seconds per gradient) x (leapfrog steps). Read from the saved fits of
## SLURM 30525 (n_c 192, reference), 30829, 30830 (43-48) and 30842 (54).
##
## Cells: per fit, summed over 4 chains. leapfrog_wu / leapfrog_s are total
## leapfrog steps (gradient evaluations) in warmup / sampling; ms_per_grad is
## total elapsed chain time divided by total leapfrog steps, in milliseconds;
## mean_lf_s is mean leapfrog steps per sampling iteration; wu_share is the
## share of chain time spent in warmup; hours is the longest chain.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc384-cost-decompose.R \
##       > _data/nc384-cost-decompose-<date>.txt

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
suppressPackageStartupMessages(library(rstan))
FITS <- c("np_bs400_nc192-seed2", "np_bs400_nc384-seed2", "np_bs15_nc384",
          "np_bs50_nc384", "np_bs100_nc384", "np_bs250_nc384", "np_bs600_nc384",
          "np_wide_total0_nc384", "np_wide_total0_nc384_t")
rows <- lapply(FITS, function(nm) {
    p <- sprintf("_data/wock-fit-%s.rds", nm)
    if (!file.exists(p)) return(NULL)
    f <- readRDS(p)
    spw <- get_sampler_params(f, inc_warmup = TRUE)
    wu <- f@stan_args[[1]]$warmup
    lf_w <- sum(sapply(spw, function(x) sum(x[seq_len(wu), "n_leapfrog__"])))
    lf_s <- sum(sapply(spw, function(x) sum(x[-seq_len(wu), "n_leapfrog__"])))
    et <- get_elapsed_time(f)
    data.frame(fit = nm,
               leapfrog_wu = lf_w, leapfrog_s = lf_s,
               ms_per_grad = round(1000 * sum(et) / (lf_w + lf_s), 1),
               mean_lf_s = round(lf_s / (4 * (f@stan_args[[1]]$iter - wu)), 1),
               wu_share = round(sum(et[, 1]) / sum(et), 2),
               hours = round(max(rowSums(et)) / 3600, 1))
})
print(do.call(rbind, rows), row.names = FALSE)
