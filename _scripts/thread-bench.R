## Time per gradient at n_c = 384 under a given thread setting. One process =
## one chain; the sbatch wrapper (_scripts/thread-bench.sh) launches several
## at once to reproduce a task's layout.
##
## Model and data are entry 43 (np_bs15_nc384), which converged; the
## gradient is evaluated at saved posterior draws from that fit, so every
## layout times identical arithmetic. log_prob is printed so that can be
## checked across layouts.
##
## Environment:
##   BENCH_KNOB     default | stan<k> | rcpp<k>   how the thread count is set:
##                  nothing; STAN_NUM_THREADS=k; RcppParallel::setThreadOptions(k)
##   BENCH_DRAW     which saved draw (one per chain, so chains differ)
##   BENCH_N        gradients timed (after 10 untimed)
##   BENCH_TAG      label for the output line
##
## Output: one line per process, key=value, to stdout.

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
knob <- Sys.getenv("BENCH_KNOB", "default")
## Set before anything can initialise the TBB pool.
if (grepl("^stan[0-9]+$", knob)) Sys.setenv(STAN_NUM_THREADS = sub("stan", "", knob))
suppressPackageStartupMessages({library(rstan); library(plasmofit)})
if (grepl("^rcpp[0-9]+$", knob))
    RcppParallel::setThreadOptions(numThreads = as.integer(sub("rcpp", "", knob)))

draw <- as.integer(Sys.getenv("BENCH_DRAW", "1"))
n    <- as.integer(Sys.getenv("BENCH_N", "200"))
tag  <- Sys.getenv("BENCH_TAG", knob)

d     <- readRDS("_data/wock-data-np_bs15_nc384.rds")
saved <- readRDS("_data/wock-fit-np_bs15_nc384.rds")
stopifnot(d$n_c == 384L, d$b_shape_fixed == 1L)

## Parameters block of archer_fit.stan for no_pool, b_shape as data, no anchor.
PARS <- c("b_off_vec", "log10_total0", "mu_logit_R", "sigma_logit_R", "eta_R",
          "mu_logit_cl", "sigma_logit_cl", "eta_cl", "z_sd_iRBC")
as_list <- function(fit, i) {
    dr <- rstan::extract(fit, pars = PARS, permuted = TRUE)
    out <- lapply(dr, function(x) { dd <- dim(x)
        if (length(dd) == 1L) x[i] else if (length(dd) == 2L) array(x[i, ], dd[2])
        else array(x[i, , ], dd[-1]) })
    for (nm in c("b_shape_free", "delta_total0", "sigma_total0"))
        out[[nm]] <- array(numeric(0), 0L)
    out
}

## A live model handle (a readRDS'd fit has no live DSO; claude/gotchas.md).
init <- as_list(saved, draw)
handle <- suppressWarnings(suppressMessages(rstan::sampling(
    plasmofit:::stanmodels[["archer_fit"]], data = d, chains = 1, iter = 1,
    algorithm = "Fixed_param", init = list(init), refresh = 0)))
u <- rstan::unconstrain_pars(handle, init)

thread_cpu <- function() {
    tt <- list.files("/proc/self/task", full.names = TRUE)
    vapply(tt, function(t) { s <- strsplit(readLines(file.path(t, "stat")), " ")[[1]]
        (as.numeric(s[14]) + as.numeric(s[15])) / 100 }, numeric(1))
}

for (k in 1:10) g <- rstan::grad_log_prob(handle, u)
c0 <- thread_cpu(); t0 <- proc.time()[["elapsed"]]
for (k in seq_len(n)) g <- rstan::grad_log_prob(handle, u)
wall <- proc.time()[["elapsed"]] - t0
c1 <- thread_cpu()
dc <- c1[names(c0)] - c0; dc[is.na(dc)] <- 0
lp <- attr(g, "log_prob")
cat(sprintf("tag=%s knob=%s draw=%d n=%d ms_per_grad=%.1f threads=%d active=%d cpu_per_wall=%.2f lp=%.10f host_cpus=%s\n",
            tag, knob, draw, n, 1000 * wall / n, length(c1), sum(dc > 0.05 * max(dc)),
            sum(dc) / wall, lp, Sys.getenv("SLURM_CPUS_PER_TASK", "NA")))
