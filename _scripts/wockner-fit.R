# Fits archer_fit.stan to the Wockner data on the cluster, one configuration
# per array task.
#
# Replaces the earlier wockner-fit.R / wock-fit-2.R / wockner-fit-cnc.R /
# cluster-test-archer-fit.R, which all hand-built the Stan data list. That job
# now belongs to plasmofit::archer_stan_data(), which also validates the input
# and fails early on things the Stan model would otherwise only reject at
# sampler launch.
#
# Add or edit entries in CONFIGS to change what gets run; the array range in
# the sbatch script below has to match its length.
#
#
# First on local machine:
#
# cd ~/GitHub/Cornell/plasmofit-test
# scp ./_scripts/wockner-fit.R ./_data/wockner-cleaned.csv lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit/
#
# Now on cluster:
#
# cd /home2/lan68/plasmofit/wock-fit
#
# ... using an interactive job:
# srun -N 1 -n 1 -c 4 --mem=8G --array=1 --time=24:00:00 --job-name="wock-fit" --pty R --vanilla
#
# ... using a non-interactive job:
#
# cat << EOF > wockner-fit.sh
# #!/bin/bash -l
#
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --array=3-7
# #SBATCH --cpus-per-task=4
# #SBATCH --mem=8G
# #SBATCH --time=2-00:00:00
# #SBATCH --job-name=wock-fit
# #SBATCH --output=wock-fit-%a.out
# #SBATCH --error=wock-fit-%a.err
# #SBATCH --mail-user=lan68@cornell.edu
# #SBATCH --mail-type=END,FAIL
#
# Rscript --vanilla wockner-fit.R
#
# EOF
#
# sbatch wockner-fit.sh
#
# Then when back on local computer:
# cd ~/GitHub/Cornell/plasmofit-test
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit/wock-fit-*.rds ./_data/
#



.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")

suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(posterior)
    library(plasmofit)
})


# ------------------------------------------------------------------------ #
# Configurations
# ------------------------------------------------------------------------ #

# Testing whether the cycle-length hierarchy earns its keep (see claude/
# CLAUDE.md "Cycle-length hierarchy: comparison in progress"). archer_fit()'s
# `model` argument (added in plasmofit's pooled-model commit ce991b1) selects
# among Stan programs that all report `log_lik` identically, so fits from
# different `model`s are directly comparable via archer_loo()/
# loo::loo_compare() -- that's the whole point of this run. "settled" data
# overrides otherwise (max_cl = 50, sd_bs_cl = 0.5); see the package's
# archer_stan_data() docs for why.
#
# `no_pool` is the baseline (current default model: R and cycle_length both
# partially pooled per grp_R/grp_cl). `pooled_cl` collapses cycle_length to a
# single value shared by every trial, leaving R hierarchical -- comparing the
# two via loo directly answers whether the cycle-length hierarchy is buying
# anything over a shared value.
#
# Both use default adapt_delta/max_treedepth. plasmofit's chat4 session found
# pooled_cl sampled badly (one chain 53% divergent) on simulated 8-series
# data at default tuning, but that was small-scale, its cause is unresolved
# (the notes and a later conversation disagree on whether the true
# cycle_length in that simulation was pooled or per-series), and real
# Wockner cycle_length only varies ~44.8-45.3 h across trials -- much
# narrower than whatever produced the simulated pathology. Check divergences/
# R-hat in the .out file once this runs; if pooled_cl struggles here too,
# rerun it with adapt_delta/max_treedepth raised rather than assuming it in
# advance.
#
# Entries 3-6 are the cycle-length prior sensitivity check, added after the
# first two runs showed the prior is doing more work than anyone intended.
# `cl_prior_center` defaults to 48 h while the data sit at ~45 h, and
# `sd_logit_cl = 1` is not weak enough for that to be ignorable: the prior
# carries weight 0.079 on no_pool's mu_logit_cl (+0.271 h) and 0.038 on
# pooled_cl's logit_cl (+0.153 h). See claude/CLAUDE.md, "cl_prior_center =
# 48 is doing real work".
#
# Two levers on the same knob: widen the prior (weight -> 0) or move its
# centre onto the data (centre - data -> 0). Widening is the cleaner
# diagnostic since it presumes no answer, so it gets both models; moving the
# centre is the confirmation that the pull is a location effect, and 42
# brackets it by pushing the opposite way.
#
# Do not widen much past 2. The prior is normal on the *logit* scale against
# hard bounds at [min_cl, max_cl]; sd_logit_cl = 3 already puts noticeable
# mass within a few hundredths of an hour of both bounds, which is how the
# max_cl = 55 boundary mode got in. The 3 run is here to find out whether
# boundary attraction comes back, not because it is wanted.
#
# Entries 1-2 are already fit and saved; submit --array=3-7 to run only the
# new ones. They are kept so the baseline stays reproducible.
CONFIGS <- list(
    no_pool        = list(model = "no_pool",   data = list()),
    pooled_cl      = list(model = "pooled_cl", data = list()),
    np_wide_prior  = list(model = "no_pool",   data = list(sd_logit_cl = 2)),
    pl_wide_prior  = list(model = "pooled_cl", data = list(sd_logit_cl = 2)),
    np_center45    = list(model = "no_pool",   data = list(cl_prior_center = 45)),
    np_center42    = list(model = "no_pool",   data = list(cl_prior_center = 42)),
    np_wider_prior = list(model = "no_pool",   data = list(sd_logit_cl = 3)),
    ## Inoculum-anchored log10_total0. `inoc` is interaction(trial, inoc_size),
    ## so inoc_size is constant within every grp_init level by construction.
    ## Predictions to read this against, from wockner-inoc-prior.R on the
    ## unanchored no_pool fit: delta_total0 ~ +0.811 (a factor of 6.47 more
    ## parasites at t = 0 than were inoculated) and R ~ 9.96.
    np_anchor      = list(model = "no_pool",   data = list(inoc_size = "inoc_size")),
    ## The control np_anchor needs. The anchored fit relaxes TWO things at
    ## once: it centres log10_total0 on the inoculum, and it replaces a tight
    ## normal(1, 0.25) with a location that is free to move (delta_total0 has
    ## sd 1). This arm relaxes only the second -- same normal centre of 1, sd
    ## widened to 1, no inoculum information at all. If it reproduces the
    ## anchored fit's elpd and its R, then the gain is the old prior being
    ## wrong and the inoculum adds nothing; if the anchored fit still wins,
    ## the inoculum carries information.
    np_wide_total0 = list(model = "no_pool",   data = list(sd_log10_total0 = 1)),
    ## b_shape gets the same treatment log10_total0 just got. lognormal(2,
    ## 0.5) has median 7.39 against a posterior of 14-19, and in simulation
    ## widening it removes 0.50 h of cycle-length bias while widening
    ## log10_total0 removes none. 10 isolates it on the original baseline; 11
    ## applies it on top of the corrected log10_total0 prior, so the two are
    ## not confounded.
    np_wide_bshape = list(model = "no_pool",   data = list(sd_log_b_shape = 1.5)),
    np_wide_both   = list(model = "no_pool",   data = list(sd_log_b_shape = 1.5,
                                                           sd_log10_total0 = 1)),
    ## The no_pool/pooled_cl hierarchy comparison was decided under the
    ## misspecified log10_total0 prior, so it has to be re-run under a
    ## defensible one before any of it is believed. Paired with 9 and 11.
    pl_wide_total0 = list(model = "pooled_cl", data = list(sd_log10_total0 = 1)),
    pl_wide_both   = list(model = "pooled_cl", data = list(sd_log_b_shape = 1.5,
                                                           sd_log10_total0 = 1)),
    ## ---- the b_shape ladder, 14-18 -------------------------------------
    ## b_shape is the shape of a symmetric Beta over cycle position, so it is
    ## a synchrony parameter and converts to a starting-stage age range (see
    ## findings.md, "What b_shape means biologically"): 14.9 -> 20.4 h on a
    ## ~45 h cycle, which is nearly no synchrony, against ~9 h assumed for
    ## controlled human infection trials in mmcm.pdf Fig. 1 (b_shape ~ 84).
    ##
    ## The panel showed b_shape is unidentified over 15-65 and its effect
    ## collapses with its value (+10 is worth 0.90 h of starting-stage sd at
    ## 15 and 0.074 h at 100), so the honest treatment is not to estimate it
    ## and not to pick one value, but to fix it across a RANGE and report the
    ## sensitivity. Each rung pins b_shape with a tight prior rather than
    ## declaring it data, which would be a package change. sd_log_b_shape =
    ## 0.05 is +-10% at 95%, which at these centres is worth <0.2 h of
    ## starting-stage sd -- pinned for every practical purpose.
    ##
    ## This is NOT the tight_sigma trap (findings.md, "Where the bias is"):
    ## that pinned a HIERARCHICAL SCALE near zero and made a funnel. This
    ## narrows a bounded scalar's own prior, which should sample no worse
    ## than the widened fits and probably better -- they carry the batch's
    ## worst divergences, 3.4-4.0%.
    ##
    ## All five carry the corrected log10_total0 prior, so the ladder is read
    ## against np_wide_total0 (task 9), not against no_pool.
    ##
    ## max_shape is raised 250 -> 400 for ALL five so the 250 rung is not
    ## sitting on the boundary -- max_cl = 55 is what reintroduced a boundary
    ## mode once already. Every rung's prior is far from 400, so the bound is
    ## non-binding throughout and the rungs stay comparable with each other.
    ## They differ from tasks 1-13 in this one respect.
    ##
    ## Prediction, recorded before the runs: relative to np_wide_total0 these
    ## should move cycle_length only ~0.1-0.3 h, because 65 -> 100 narrows the
    ## age range 10.1 -> 8.2 h where 15 -> 65 narrowed it 20.4 -> 10.1 h. A
    ## much larger move means this reasoning is wrong.
    np_bs50        = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(50),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 400)),
    np_bs84        = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(84),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 400)),
    np_bs100       = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(100),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 400)),
    np_bs150       = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(150),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 400)),
    np_bs250       = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(250),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 400)),
    ## ---- extending the ladder, 19-21 -----------------------------------
    ## elpd was still rising monotonically at the 250 rung (+36.2 against a
    ## free b_shape, and +2.26 against the 100 rung), so 250 is the top of
    ## what was run rather than an optimum. These bracket it.
    ##
    ## 400 and 600 need max_shape above them; 1000 is used for both. That
    ## makes them differ from rungs 14-18 in two things at once, since
    ## max_shape also sets the scale of Stan's bounded transform and so can
    ## move divergences and step size on its own -- and divergences rising
    ## with the pin is one of the ladder's findings. np_bs250_ms1000 is the
    ## control that separates the two: it is rung 18 with max_shape alone
    ## changed, so the difference between it and np_bs250 is what max_shape
    ## is worth, and anything beyond that in 400/600 is the pin.
    ##
    ## An age range of 5.2 h at 250 falls to 4.1 h at 400 and 3.4 h at 600,
    ## against 0.47 h per stage cell -- still several cells wide, so this is
    ## not yet asking the model to resolve below its own discretisation.
    np_bs400       = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(400),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 1000)),
    np_bs600       = list(model = "no_pool",   data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(600),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 1000)),
    np_bs250_ms1000 = list(model = "no_pool",  data = list(sd_log10_total0 = 1,
                                                           mean_log_b_shape = log(250),
                                                           sd_log_b_shape = 0.05,
                                                           max_shape = 1000))
)

# log_lik is needed for loo/waic but roughly triples the size of a stored fit.
# Needed here since the whole point of this run is a loo comparison between
# models -- this is also the first saved fit with it (see CLAUDE.md).
CALC_LOG_LIK <- TRUE

# Fixed, because this posterior is multimodal: without it, two runs differ by
# which mode each chain initializes into and not just by Monte Carlo error.
SEED <- 538065874

N_CHAINS <- 4L
WARMUP <- 700L
ITER <- 1000L + WARMUP


n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")
if (curr_idx > length(CONFIGS)) {
    stop("SLURM_ARRAY_TASK_ID ", curr_idx, " exceeds the ", length(CONFIGS),
         " entries in CONFIGS")
}

cfg_name <- names(CONFIGS)[curr_idx]
cfg <- CONFIGS[[curr_idx]]

# Within-chain threading measured at roughly 1.0x on these data: there are only
# ~14 distinct trajectory combinations to split over, so the scheduler has
# little to work with. Chains are the better use of cores here. Raise this only
# for a dataset with many more combinations, and measure before believing it.
threads_per_chain <- 1L
options(mc.cores = min(N_CHAINS, max(1L, n_threads %/% threads_per_chain)))


# ------------------------------------------------------------------------ #
# Data
# ------------------------------------------------------------------------ #

paras_df <- read_csv("wockner-cleaned.csv", col_types = "cccdcdd") |>
    # Group for inoculations:
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |> paste()) |>
    # Group for observation error:
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste())

d <- do.call(archer_stan_data,
             c(list(data = paras_df,
                    series = "id", time = "time", abundance = "para",
                    grp_init = "inoc", grp_R = "trial", grp_cl = "trial",
                    grp_sd = "obs_error",
                    calc_log_lik = as.integer(CALC_LOG_LIK)),
               cfg$data))

cat("config:", cfg_name, "| model:", cfg$model, "\n")
if (length(cfg$data) > 0) {
    cat("  data overrides:", paste(names(cfg$data), unlist(cfg$data),
                                   sep = " = ", collapse = ", "), "\n")
}
if (length(cfg$control) > 0) {
    cat("  control:", paste(names(cfg$control), unlist(cfg$control),
                            sep = " = ", collapse = ", "), "\n")
}
cat(sprintf("  %d series, %d observations, n_c = %d, dt_full = %g\n",
            d$n_ts, d$n_total_obs, d$n_c, d$dt_full))
cat(sprintf("  groups init/R/cl/sd: %d/%d/%d/%d\n",
            d$n_grp_init, d$n_grp_R, d$n_grp_cl, d$n_grp_sd))
cat(sprintf("  cycle_length in (%g, %g), prior centred at %.2f h, sd_bs_cl %g\n",
            d$min_cl, d$max_cl,
            d$min_cl + (d$max_cl - d$min_cl) * inv_logit(d$mean_logit_cl),
            d$sd_bs_cl))


# ------------------------------------------------------------------------ #
# Fit
# ------------------------------------------------------------------------ #

f <- archer_fit(d, model = cfg$model, chains = N_CHAINS, iter = ITER,
                warmup = WARMUP, seed = SEED,
                threads_per_chain = threads_per_chain, control = cfg$control)

write_rds(f, sprintf("wock-fit-%s.rds", cfg_name))
write_rds(d, sprintf("wock-data-%s.rds", cfg_name))

cat("\nFinished", cfg_name, "\n")
print(rstan::get_elapsed_time(f))
cat("\n-------------------------------------\n")
cat(cfg_name,
    "| leapfrog", round(mean(rstan::get_num_leapfrog_per_iteration(f)), 1),
    "| div", sum(rstan::get_divergent_iterations(f)), "\n")


## loo objects, saved rather than just printed: loo_compare() needs the
## objects themselves, and recomputing them needs log_lik, which we may not
## keep. Done first, right after f/d are saved, so the actual point of this
## run (the loo comparison) survives even if a diagnostic below turns out
## fragile for a model it wasn't written for.
if (CALC_LOG_LIK && requireNamespace("loo", quietly = TRUE)) {
    loos <- list(observation = archer_loo(f, d),
                 trial       = archer_loo(f, d, by = "grp_cl"))
    write_rds(loos, sprintf("wock-fit-LOO-%s.rds", cfg_name))
    for (nm in names(loos)) {
        k <- loo::pareto_k_values(loos[[nm]])
        cat(sprintf("\nloo (%s): elpd %.1f (se %.1f), %d units, %d with k > 0.7\n",
                    nm, loos[[nm]]$estimates["elpd_loo", "Estimate"],
                    loos[[nm]]$estimates["elpd_loo", "SE"], length(k), sum(k > 0.7)))
    }
    cat("\nNote: leaving out a whole trial perturbs the posterior much more than\n",
        "leaving out one observation, so high Pareto k here is expected and means\n",
        "the approximation, not the model, is failing. A trustworthy answer at\n",
        "that grouping needs K-fold refitting.\n", sep = "")
}


# ------------------------------------------------------------------------ #
# Diagnostics that need a live DSO, so they have to happen here rather than
# after the fit is written out and reloaded elsewhere.
#
# Which R/cycle_length parameters exist depends on `cfg$model`: pooled_R
# replaces mu_logit_R/sigma_logit_R/eta_R with a single logit_R, and
# pooled_cl does the same for cycle_length (mu_logit_cl/sigma_logit_cl/
# eta_cl -> logit_cl). Generated quantities (R, cycle_length, b_offset, ...)
# are identically shaped across all four models, so only the code that
# touches parameters directly needs to branch on the model.
# ------------------------------------------------------------------------ #

r_pars_for  <- function(model) {
    if (model %in% c("pooled_R", "pooled_both")) "logit_R"
    else c("mu_logit_R", "sigma_logit_R", "eta_R")
}
cl_pars_for <- function(model) {
    if (model %in% c("pooled_cl", "pooled_both")) "logit_cl"
    else c("mu_logit_cl", "sigma_logit_cl", "eta_cl")
}
## The inoculum anchor adds two free parameters, declared array[total0_anchor]
## and so zero-sized when it is off. unconstrain_pars() needs EVERY free
## parameter, so omitting them would fail on an anchored fit -- after the
## sampling is already paid for.
anchor_pars_for <- function(data) {
    if (isTRUE(as.integer(data$total0_anchor) == 1L))
        c("delta_total0", "sigma_total0") else character(0)
}
## Since b_shape became optionally-data, the PARAMETER is `b_shape_free` and
## `b_shape` is a transformed parameter assembled from it or from the data.
## unconstrain_pars() wants the parameter, and passing it a transformed one
## errors. Models compiled before that change still have `b_shape` as the
## parameter, so pick whichever this fit actually has -- fit@model_pars lists
## transformed parameters too, so the presence of `b_shape` cannot settle it
## and `b_shape_free` has to be the thing tested for.
bshape_par_for <- function(fit) {
    if ("b_shape_free" %in% fit@model_pars) "b_shape_free" else "b_shape"
}

par_block_for <- function(model, fit, data = d) {
    c(bshape_par_for(fit), "b_off_vec", "log10_total0",
      r_pars_for(model), cl_pars_for(model), "z_sd_iRBC",
      anchor_pars_for(data))
}

## one draw as a named list shaped the way unconstrain_pars expects.
##
## unconstrain_pars() needs an entry for EVERY parameter the model declares,
## including ones declared array[0] and therefore zero-sized. Since the
## inoculum anchor was added, archer_fit.stan always declares delta_total0 and
## sigma_total0, so a fit with the anchor OFF still needs them present at
## length 0 -- otherwise "variable does not exist". Fits compiled before that
## change do not declare them at all, so they are added only when the fit says
## it has them.
draw_as_list <- function(fit, pars, i = 1L) {
    dr <- rstan::extract(fit, pars = pars, permuted = TRUE)
    out <- lapply(dr, function(x) {
        d <- dim(x)
        if (length(d) == 1L) x[i]
        ## array(), not x[i, ]: a length-one container (array[1] real, as the
        ## anchor's delta_total0 and sigma_total0 are declared) would otherwise
        ## come back as a bare scalar with no dim, and unconstrain_pars() reads
        ## dims from the object -- "dims declared=(1); dims found=()".
        else if (length(d) == 2L) array(x[i, ], d[2])
        else array(x[i, , ], d[-1])
    })
    ## Zero-sized parameters contribute no draws, so rstan::extract returns
    ## nothing for them -- but unconstrain_pars() still needs an entry or it
    ## fails with "variable does not exist". b_shape_free is zero-sized
    ## whenever b_shape was supplied as data, exactly as delta_total0 and
    ## sigma_total0 are whenever the anchor is off.
    for (nm in intersect(c("delta_total0", "sigma_total0", "b_shape_free"),
                         fit@model_pars)) {
        if (!nm %in% names(out)) out[[nm]] <- array(numeric(0), 0L)
    }
    out
}

## seconds per gradient evaluation
grad_time <- function(fit, model, n = 200L) {
    u <- rstan::unconstrain_pars(fit, draw_as_list(fit, par_block_for(model, fit), 1L))
    system.time(for (k in seq_len(n)) rstan::grad_log_prob(fit, u))[["elapsed"]] / n
}

## condition number of the posterior correlation matrix, unconstrained scale
fit_cond <- function(fit, max_shape, model) {
    lg <- function(p) log(p / (1 - p))
    pars <- par_block_for(model, fit)
    dr <- rstan::extract(fit, pars = pars, permuted = TRUE)
    ## b_shape contributes a block only when it was ESTIMATED. Supplied as
    ## data it is constant, so it has no unconstrained coordinate and a
    ## constant column would make cor(U) undefined (zero variance -> NaN)
    ## rather than merely uninformative.
    bs_nm <- bshape_par_for(fit)
    bs <- if (!is.null(dr[[bs_nm]]) && length(dr[[bs_nm]]) > 0) {
        lg((dr[[bs_nm]] - 2) / (max_shape - 2))
    } else NULL
    ## unit_vector[2] is rank-deficient in storage; use the angle instead
    ang <- atan2(dr$b_off_vec[, , 2], dr$b_off_vec[, , 1])
    colnames(ang) <- paste0("b_ang[", seq_len(ncol(ang)), "]")
    ## logit_R/logit_cl are already unconstrained scalars, so they need no
    ## further transform, unlike mu_logit_*/log(sigma_logit_*)/eta_*
    r_block <- if ("logit_R" %in% pars) {
        matrix(dr$logit_R, ncol = 1, dimnames = list(NULL, "logit_R"))
    } else {
        cbind(mu_logit_R = dr$mu_logit_R, sigma_logit_R = log(dr$sigma_logit_R), dr$eta_R)
    }
    cl_block <- if ("logit_cl" %in% pars) {
        matrix(dr$logit_cl, ncol = 1, dimnames = list(NULL, "logit_cl"))
    } else {
        cbind(mu_logit_cl = dr$mu_logit_cl, sigma_logit_cl = log(dr$sigma_logit_cl), dr$eta_cl)
    }
    ## cbind() drops a NULL, so a fixed b_shape simply leaves the block out.
    U <- cbind(bs, ang, dr$log10_total0, r_block, cl_block, dr$z_sd_iRBC)
    ev <- eigen(cor(U), only.values = TRUE)$values
    sqrt(max(ev) / min(ev))
}

## do divergences concentrate where a scale parameter is small? that is the
## signature of a funnel rather than of generic trouble. A pooled model has
## no sigma for the parameter(s) it pools, so there may be nothing to check.
div_by_sigma <- function(fit, model) {
    div <- as.logical(rstan::get_divergent_iterations(fit))
    if (!any(div)) return(NULL)
    sigma_pars <- c(if (!(model %in% c("pooled_R", "pooled_both"))) "sigma_logit_R",
                    if (!(model %in% c("pooled_cl", "pooled_both"))) "sigma_logit_cl")
    if (length(sigma_pars) == 0) return(NULL)
    s <- as.matrix(fit, pars = sigma_pars)
    rbind(divergent = colMeans(s[div, , drop = FALSE]),
          ok        = colMeans(s[!div, , drop = FALSE]))
}

summarize_fit <- function(fit, data, label, model,
                          par_names = c("b_shape", "b_offset", "R",
                                        "log10_total0", "cycle_length",
                                        "sd_iRBC")) {
    post <- rstan::extract(fit, permuted = FALSE)
    diags <- do.call(rbind, lapply(par_names, function(p) {
        sims <- rstan::extract(fit, p, permuted = FALSE)
        data.frame(par      = dimnames(sims)$parameters,
                   rhat     = apply(sims, 3, posterior::rhat),
                   ess_bulk = apply(sims, 3, posterior::ess_bulk),
                   ess_tail = apply(sims, 3, posterior::ess_tail),
                   row.names = NULL)
    }))
    ## drop = FALSE: with a single chain this otherwise collapses to a vector
    ## and ncol() comes back NULL
    lp <- post[, , "lp__", drop = FALSE]
    list(
        label        = label,
        model        = model,
        center_R     = data$center_R,
        center_cl    = data$center_cl,
        max_cl       = data$max_cl,
        sd_bs_cl     = data$sd_bs_cl,
        fit          = fit,
        lp_by_chain  = apply(lp, 2, mean),
        leapfrog     = mean(rstan::get_num_leapfrog_per_iteration(fit)),
        n_divergent  = sum(rstan::get_divergent_iterations(fit)),
        div_sigma    = div_by_sigma(fit, model),
        stepsize     = sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                              function(x) mean(x[, "stepsize__"])),
        elapsed      = rstan::get_elapsed_time(fit),
        sec_per_grad = grad_time(fit, model),
        cond         = fit_cond(fit, max_shape = data$max_shape, model = model),
        diagnostics  = diags
    )
}

res <- summarize_fit(f, d, label = cfg_name, model = cfg$model)
write_rds(res, sprintf("wock-fit-RES-%s.rds", cfg_name))
