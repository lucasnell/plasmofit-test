# Is the cycle-length bias a MARGINALISATION effect?
#
# The sharpest fact in findings.md is that the likelihood is not biased and
# the posterior is: the per-trial MLE recovers the truth to -0.195 h, the
# pooled MLE to +0.544 h, and the fitted posterior mean sits at +1.56 to
# +1.71 h. Widening the nuisance priors buys only ~0.43 h of that, so the
# standing attribution -- "nuisance priors dragging cycle_length along the
# ridge" -- does not cover the gap.
#
# Mean, median and marginal mode cannot distinguish the possibilities,
# because all three are MARGINAL summaries: each integrates over the other
# ~200 parameters alike, so a bias created by that integration appears in all
# of them equally (measured: they span 0.05 h on a 1.8 h bias). A JOINT
# summary does not integrate at all, so comparing one against the posterior
# mean isolates the integration.
#
# rstan::optimizing maximises the joint posterior on the constrained scale,
# which makes a three-point ladder on the SAME simulated data:
#
#   pooled MLE      likelihood only, sd_iRBC fixed at truth     +0.544 h
#   penalised max   likelihood + every prior the fit uses       <- this script
#   posterior mean  the same thing, marginalised                +1.56 to +1.71 h
#
# Near the MLE  -> the gap is marginalisation, and the question becomes which
#                  summary to report rather than what is broken.
# Near the mean -> the gap is in the priors jointly, even though widening the
#                  two nuisance ones individually did not move it, and the
#                  cycle-length prior and the hierarchy are then the suspects.
#
# MULTIPLE RANDOM STARTS, because this posterior is multimodal -- it is why
# wockner-fit.R fixes its seed. A single optimisation would report whichever
# mode it happened to fall into. The best start by lp__ is the estimate; the
# spread across starts is reported beside it, and if the starts disagree
# materially the MAP is not a well-defined target and that is the finding.
#
#   SCHEDSIM_TASK=<n> picks the replicate: 1, 2, 3 = default-rep1..3.
#   MAP_N_INIT=<n> sets the number of random starts (default 12).
#
#     SCHEDSIM_TASK=2 srun -N 1 -n 1 -c 4 --mem=32G Rscript --vanilla \
#         _scripts/schedsim-map-check.R

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

TASK <- as.integer(Sys.getenv("SCHEDSIM_TASK", "2"))
N_INIT <- as.integer(Sys.getenv("MAP_N_INIT", "12"))

## Rebuild the simulated data exactly as the fit saw it, by running
## wockner-schedule-sim.R's own prologue -- deterministic from the recorded
## noise seed, so this is the same dataset and not a re-simulation of it.
## Sourced rather than copied: a copy would drift from the script that made
## the fits, which is the whole premise of the comparison.
Sys.setenv(SLURM_ARRAY_TASK_ID = as.character(TASK))
src <- readLines("_scripts/wockner-schedule-sim.R")
cut <- grep("^# Fit$", src) - 2L
eval(parse(text = paste(src[seq_len(cut)], collapse = "\n")))

suppressPackageStartupMessages({ library(tidyverse); library(rstan) })

cat("\n=== MAP check:", cfg_name, "===\n")
cat("true cycle_length:", truth$cycle_length, "h | random starts:", N_INIT, "\n\n")

sm <- plasmofit:::stanmodels[[if (arm$model == "no_pool") "archer_fit"
                              else paste0("archer_fit_", arm$model)]]

runs <- map(seq_len(N_INIT), \(i) {
    t0 <- Sys.time()
    o <- tryCatch(
        rstan::optimizing(sm, data = d_sim, seed = 1000L + i,
                          init = "random", as_vector = FALSE,
                          hessian = FALSE, verbose = FALSE),
        error = function(e) list(return_code = -1L, message = conditionMessage(e)))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (!is.null(o$return_code) && o$return_code == 0L) {
        cl <- o$par$cycle_length
        cat(sprintf("  start %2d: lp__ %12.2f  mean cycle_length %.4f  (%.1f min)\n",
                    i, o$value, mean(cl), el))
        ## b_offset is the phase. An oscillator's period and phase trade off
        ## -- a period that is slightly wrong can be compensated by shifting
        ## the phase -- so if the modes lie on a (period, phase) comb rather
        ## than scattering, the multimodality is that trade-off and not
        ## generic roughness. Recorded per start so the question is
        ## answerable; the first version of this script saved only
        ## cycle_length and could not tell the two apart.
        tibble(start = i, ok = TRUE, lp = o$value, mean_cl = mean(cl),
               min_cl = min(cl), max_cl = max(cl),
               mean_boff = mean(o$par$b_offset),
               boff1 = o$par$b_offset[1],
               mean_bshape = mean(o$par$b_shape),
               mean_R = mean(o$par$R),
               mean_sd = mean(o$par$sd_iRBC),
               norm_boff = mean(sqrt(rowSums(o$par$b_off_vec^2))))
    } else {
        cat(sprintf("  start %2d: FAILED (%s)\n", i,
                    o$message %||% paste("return_code", o$return_code)))
        tibble(start = i, ok = FALSE, lp = NA_real_, mean_cl = NA_real_,
               min_cl = NA_real_, max_cl = NA_real_, mean_boff = NA_real_,
               boff1 = NA_real_, mean_bshape = NA_real_, mean_R = NA_real_,
               mean_sd = NA_real_, norm_boff = NA_real_)
    }
}) |> list_rbind()

ok <- runs |> filter(ok)
cat("\nCells: one row per random start that converged. `lp` is the joint\n")
cat("log posterior at the optimum, so the largest is the best mode found;\n")
cat("`mean_cl` is cycle_length averaged over the 13 trials at that optimum,\n")
cat("in hours. Starts are independent; spread across them measures how\n")
cat("multimodal this surface is, not Monte Carlo error.\n\n")
runs |> filter(ok) |> slice_max(lp, n = 15) |>
    mutate(across(where(is.numeric), \(x) round(x, 4))) |> print(n = Inf, width = Inf)
cat("\n(top 15 by lp; all starts are in the saved rds)\n")

if (nrow(ok) > 0) {
    best <- ok |> slice_max(lp, n = 1)
    cat(sprintf("\nbest mode: lp__ %.2f, mean cycle_length %.4f h, bias %+.4f h\n",
                best$lp, best$mean_cl, best$mean_cl - truth$cycle_length))
    cat(sprintf("across %d converged starts: lp range %.2f, cycle_length range %.4f h\n",
                nrow(ok), diff(range(ok$lp)), diff(range(ok$mean_cl))))
    ## how many starts are within a trivial lp of the best, i.e. the same mode
    near <- sum(ok$lp > best$lp - 0.01)
    cat(sprintf("starts reaching the best mode (lp within 0.01): %d of %d\n",
                near, nrow(ok)))
}

write_rds(list(config = cfg_name, truth_cl = truth$cycle_length,
               n_init = N_INIT, runs = runs),
          sprintf("_data/wock-map-check-%s.rds", cfg_name))
cat("\nwrote _data/wock-map-check-", cfg_name, ".rds\n", sep = "")
