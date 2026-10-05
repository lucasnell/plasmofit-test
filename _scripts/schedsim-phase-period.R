# Is the cycle-length bias a phase-period trade-off?
#
# Thread 2 is out of named suspects. This tests an unnamed one, suggested by
# two numbers already in findings.md's recovery table that have never been
# read together: `cycle_length` recovers +1.9 h and `b_offset` recovers +45%,
# i.e. +0.115 cycles. Those are not independent errors. A longer period loses
# phase over the observation window -- 216 h is 4.80 cycles at the truth and
# 4.61 at the estimate, a loss of 0.193 cycles -- and a larger initial offset
# gains it back. The observed b_offset error recovers 60% of exactly that.
#
# If the data constrain the PHASE AT OBSERVATION TIMES and only weakly
# separate the period from the initial offset, then cycle_length's bias is
# not an error the data could have corrected: it is one coordinate of a ridge
# the likelihood is nearly flat along, and the posterior mean lands wherever
# the ridge geometry and the prior put it.
#
# The test: compare how well phi(t) = frac(b_offset + t / cycle_length) is
# recovered against how well its two components are. phi is the natural
# reading of this parameterisation -- beta_starts() shifts the initial stage
# distribution by b_offset and the trajectory advances with cycle_length --
# and it is stated here as a reading, not derived from the Stan code.
#
# Errors are CIRCULAR: phase lives on [0, 1), so the difference is wrapped to
# [-0.5, 0.5). Reporting a raw difference would make a 0.99 vs 0.01 error
# look like 0.98 instead of 0.02.
#
#   SCHEDSIM_TASK=<n> picks the replicate (1, 2, 3 = default-rep1..3).
#     SCHEDSIM_TASK=1 srun -N 1 -n 1 -c 4 --mem=32G Rscript --vanilla \
#         _scripts/schedsim-phase-period.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
TASK <- as.integer(Sys.getenv("SCHEDSIM_TASK", "1"))
Sys.setenv(SLURM_ARRAY_TASK_ID = as.character(TASK))
src <- readLines("_scripts/wockner-schedule-sim.R")
cut <- grep("^# Fit$", src) - 2L
eval(parse(text = paste(src[seq_len(cut)], collapse = "\n")))

suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

wrap <- function(x) ((x + 0.5) %% 1) - 0.5      # to [-0.5, 0.5)

## Phase is CIRCULAR and its posterior here is wide -- a quarter of a cycle
## on some series -- so an arithmetic mean of wrapped differences is pulled
## toward zero and is not the centre of anything. The first version of this
## script used one and reported phase errors roughly a third too small.
## Circular mean and circular sd, both in cycles.
cmean <- function(x) atan2(mean(sin(2 * pi * x)), mean(cos(2 * pi * x))) / (2 * pi)
csd <- function(x) {
    R <- sqrt(mean(sin(2 * pi * x))^2 + mean(cos(2 * pi * x))^2)
    sqrt(-2 * log(max(R, 1e-12))) / (2 * pi)
}

f <- read_rds(sprintf("_data/wock-schedsim-fit-%s.rds", cfg_name))
CL <- as.matrix(f, pars = "cycle_length")       # draws x n_grp_cl
BO <- as.matrix(f, pars = "b_offset")           # draws x n_grp_init
rm(f); invisible(gc())

ends <- cumsum(d$n_obs); starts <- ends - d$n_obs + 1L
t_first <- d$ts[starts]; t_last <- d$ts[ends]

cat("\n=== replicate:", cfg_name, "===\n")
cat("true cycle_length", round(truth$cycle_length, 4),
    "| mean true b_offset", round(mean(truth$b_offset), 4), "\n")
cat("observation span: first", min(t_first), "h, last", max(t_last), "h =",
    round(max(t_last) / truth$cycle_length, 2), "cycles at the truth\n\n")

phi <- function(bo, cl, t) (bo + t / cl) %% 1

## per series, averaged over draws then over series
per_series <- map(seq_len(d$n_ts), \(i) {
    gi <- d$grp_init[i]; gc <- d$grp_cl[i]
    cl_d <- CL[, gc]; bo_d <- BO[, gi]
    tibble(series = i,
           err_cl    = mean(cl_d) - truth$cycle_length,
           err_bo    = wrap(mean(bo_d) - truth$b_offset[gi]),
           err_first = cmean(phi(bo_d, cl_d, t_first[i]) -
                             phi(truth$b_offset[gi], truth$cycle_length, t_first[i])),
           err_last  = cmean(phi(bo_d, cl_d, t_last[i]) -
                             phi(truth$b_offset[gi], truth$cycle_length, t_last[i])),
           sd_last   = csd(phi(bo_d, cl_d, t_last[i])),
           cor_cl_bo = cor(cl_d, bo_d))
}) |> list_rbind()

cat("Cells: averaged over the", d$n_ts, "series. `err_cl` is the posterior
mean cycle_length minus the truth, in HOURS. The other errors are in CYCLES,
wrapped to [-0.5, 0.5): `err_bo` for the initial offset, `err_first` and
`err_last` for the phase at each series' first and last observation. `cl in
cycles` converts err_cl to the phase it costs over the mean observation span,
which is the quantity b_offset would have to cancel. `sd_last` is the
posterior CIRCULAR sd of the phase at the last observation, in cycles -- how
tightly the data pin the phase at all. `cor` is the within-draw posterior
correlation between a series' cycle_length and its b_offset. Phase summaries
are circular means; an arithmetic mean on a quantity this diffuse is not a
centre.\n\n")

span <- mean(t_last)
summ <- per_series |>
    summarise(err_cl_h = mean(err_cl),
              cl_in_cycles = mean(-err_cl / truth$cycle_length *
                                      (span / truth$cycle_length)),
              err_bo = mean(err_bo),
              err_first = cmean(err_first), err_last = cmean(err_last),
              sd_last = mean(sd_last),
              cor = mean(cor_cl_bo))
summ |> mutate(across(everything(), \(x) round(x, 4))) |> print(width = Inf)

cat(sprintf("\n  |phase error at the LAST observation| = %.4f cycles = %.2f h\n",
            abs(summ$err_last), abs(summ$err_last) * truth$cycle_length))
cat(sprintf("  |cycle_length error|                 = %.4f cycles = %.2f h\n",
            abs(summ$cl_in_cycles), abs(summ$err_cl_h)))
cat(sprintf("  ratio: the phase at the last observation is recovered %.1fx\n",
            abs(summ$cl_in_cycles) / max(abs(summ$err_last), 1e-9)))
cat("  better than the period that supposedly produces it.\n")
