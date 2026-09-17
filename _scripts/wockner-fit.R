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
# #SBATCH --array=1-3
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

# Each entry is a set of overrides passed to archer_stan_data(). Anything not
# named here takes the package default, which is the settled configuration:
# max_cl = 50 and sd_bs_cl = 0.5.
#
# Those two were settled empirically. max_cl = 55 admits a second mode at the
# bound (~54 h) that fits ~16 log units worse and wrecks R-hat; sd_bs_cl = 0.1
# is tight enough that the prior, not the data, sets the between-trial spread,
# and it pushes sigma_logit_cl into the funnel neck. The two older variants are
# kept below so those claims stay reproducible rather than being folk memory.
CONFIGS <- list(
    settled      = list(),
    wide_max_cl  = list(max_cl = 55),   # reintroduces the boundary mode
    tight_sd_bs  = list(sd_bs_cl = 0.1) # prior-dominated between-trial spread
)

# log_lik is needed for loo/waic but roughly triples the size of a stored fit,
# so keep it on only where a model comparison is actually planned.
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
               cfg))

cat("config:", cfg_name, "\n")
if (length(cfg) > 0) {
    cat("  overrides:", paste(names(cfg), unlist(cfg), sep = " = ",
                              collapse = ", "), "\n")
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

f <- archer_fit(d, chains = N_CHAINS, iter = ITER, warmup = WARMUP,
                seed = SEED, threads_per_chain = threads_per_chain)

write_rds(f, sprintf("wock-fit-%s.rds", cfg_name))
write_rds(d, sprintf("wock-data-%s.rds", cfg_name))

cat("\nFinished", cfg_name, "\n")
print(rstan::get_elapsed_time(f))
cat("\n-------------------------------------\n")
cat(cfg_name,
    "| leapfrog", round(mean(rstan::get_num_leapfrog_per_iteration(f)), 1),
    "| div", sum(rstan::get_divergent_iterations(f)), "\n")


# ------------------------------------------------------------------------ #
# Diagnostics that need a live DSO, so they have to happen here rather than
# after the fit is written out and reloaded elsewhere.
# ------------------------------------------------------------------------ #

PAR_BLOCK <- c("b_shape", "b_off_vec", "log10_total0",
               "mu_logit_R", "sigma_logit_R", "eta_R",
               "mu_logit_cl", "sigma_logit_cl", "eta_cl",
               "z_sd_iRBC")

## one draw as a named list shaped the way unconstrain_pars expects
draw_as_list <- function(fit, i = 1L, pars = PAR_BLOCK) {
    dr <- rstan::extract(fit, pars = pars, permuted = TRUE)
    lapply(dr, function(x) {
        d <- dim(x)
        if (length(d) == 1L) x[i]
        else if (length(d) == 2L) x[i, ]
        else array(x[i, , ], d[-1])
    })
}

## seconds per gradient evaluation
grad_time <- function(fit, n = 200L) {
    u <- rstan::unconstrain_pars(fit, draw_as_list(fit, 1L))
    system.time(for (k in seq_len(n)) rstan::grad_log_prob(fit, u))[["elapsed"]] / n
}

## condition number of the posterior correlation matrix, unconstrained scale
fit_cond <- function(fit, max_shape) {
    lg <- function(p) log(p / (1 - p))
    dr <- rstan::extract(fit, pars = PAR_BLOCK, permuted = TRUE)
    bs <- lg((dr$b_shape - 2) / (max_shape - 2))
    ## unit_vector[2] is rank-deficient in storage; use the angle instead
    ang <- atan2(dr$b_off_vec[, , 2], dr$b_off_vec[, , 1])
    colnames(ang) <- paste0("b_ang[", seq_len(ncol(ang)), "]")
    U <- cbind(bs, ang, dr$log10_total0,
               mu_logit_R     = dr$mu_logit_R,
               sigma_logit_R  = log(dr$sigma_logit_R),
               dr$eta_R,
               mu_logit_cl    = dr$mu_logit_cl,
               sigma_logit_cl = log(dr$sigma_logit_cl),
               dr$eta_cl,
               dr$z_sd_iRBC)
    ev <- eigen(cor(U), only.values = TRUE)$values
    sqrt(max(ev) / min(ev))
}

## do divergences concentrate where a scale parameter is small? that is the
## signature of a funnel rather than of generic trouble
div_by_sigma <- function(fit) {
    div <- as.logical(rstan::get_divergent_iterations(fit))
    if (!any(div)) return(NULL)
    s <- as.matrix(fit, pars = c("sigma_logit_R", "sigma_logit_cl"))
    rbind(divergent = colMeans(s[div, , drop = FALSE]),
          ok        = colMeans(s[!div, , drop = FALSE]))
}

summarize_fit <- function(fit, data, label,
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
        center_R     = data$center_R,
        center_cl    = data$center_cl,
        max_cl       = data$max_cl,
        sd_bs_cl     = data$sd_bs_cl,
        fit          = fit,
        lp_by_chain  = apply(lp, 2, mean),
        leapfrog     = mean(rstan::get_num_leapfrog_per_iteration(fit)),
        n_divergent  = sum(rstan::get_divergent_iterations(fit)),
        div_sigma    = div_by_sigma(fit),
        stepsize     = sapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                              function(x) mean(x[, "stepsize__"])),
        elapsed      = rstan::get_elapsed_time(fit),
        sec_per_grad = grad_time(fit),
        cond         = fit_cond(fit, max_shape = data$max_shape),
        diagnostics  = diags
    )
}

res <- summarize_fit(f, d, label = cfg_name)
write_rds(res, sprintf("wock-fit-RES-%s.rds", cfg_name))


## loo objects, saved rather than just printed: loo_compare() needs the objects
## themselves, and recomputing them needs log_lik, which we may not keep.
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
