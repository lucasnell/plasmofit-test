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
# srun -N 1 -n 1 -c 4 --mem=40G --array=2 --time=24:00:00 --job-name="wock-fit" --pty R --vanilla
#
# ... using a non-interactive job:
#
# cat << EOF > wockner-fit.sh
# #!/bin/bash -l
#
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --array=1-4
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
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")


suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(posterior)
    library(plasmofit)
})

.n_chains <- 4L

rstan_options(threads_per_chain = max(1L, n_threads %/% .n_chains))
options(mc.cores = max(1L, n_threads %/% rstan_options("threads_per_chain")))

stopifnot((options()[["mc.cores"]] * rstan_options("threads_per_chain")) <= n_threads)


# w <- c(400L, 700L, 1000L)[curr_idx]
w <- 700L


# Which change to make (or no change if "none"):
change <- list(none = NA_real_,
               sd_bs_cl = 0.5,
               max_cl = 50,
               center_cl = 0)[curr_idx]
stopifnot(length(change) == 1L)


paras_df <- read_csv("wockner-cleaned.csv", col_types = "cccdcdd") |>
    # Group for inoculations:
    mutate(inoc = interaction(trial, inoc_size, drop = TRUE) |>
               paste()) |>
    # Group for observation error:
    mutate(obs_error = interaction(trial, cohort, drop = TRUE) |> paste()) |>
    # Order these factors by order of appearance:
    mutate(across(all_of(c("id", "trial", "inoc", "obs_error")),
                  \(x) factor(x, levels = unique(x))))


summ_by_ts_paras_df <- paras_df |>
    group_by(id) |>
    summarize(inoc = inoc[1],
              trial = trial[1],
              obs_error = obs_error[1]) |>
    mutate(across(inoc:obs_error, as.integer))


# greatest common denominator (assumes all x and y are integers)
gcd <- function(x, y) ifelse(y == 0, x, gcd(y, x %% y))


stan_env <- new.env()

with(stan_env, {

    # Numbers of groups:
    n_ts <- length(levels(paras_df$id))
    n_grp_init <- length(levels(paras_df$inoc))
    n_grp_R <- length(levels(paras_df$trial))
    n_grp_cl <- length(levels(paras_df$trial))
    n_grp_sd <- length(levels(paras_df$obs_error))

    # Data:
    ts <- paras_df$time
    y <- paras_df$para

    # Other integers:
    n_obs <- paras_df |>
        group_by(id) |>
        summarize(n_obs = n()) |>
        getElement("n_obs")
    n_total_obs <- sum(n_obs)
    mu <- 0
    n_c <- 96L
    dt_full <- Reduce(gcd, paras_df$time)
    run_check <- 0L
    calc_log_lik <- 0L   # 1L to output pointwise log_lik for loo/waic
    center_R <- 1L
    center_cl <- 1L
    grainsize <- 1L

    # Indices for grouping time series with the same for each:
    grp_init <- summ_by_ts_paras_df$inoc
    grp_R <- summ_by_ts_paras_df$trial
    grp_cl <- summ_by_ts_paras_df$trial
    grp_sd <- summ_by_ts_paras_df$obs_error

    # Parameter limitations:
    max_shape <- 250
    max_R <- 50
    min_cl <- 35
    max_cl <- 50

    # Hyperparameters:
    mean_log_b_shape <- rep(2, n_grp_init)
    sd_log_b_shape <- rep(0.5, n_grp_init)
    mean_log10_total0 <- rep(1, n_grp_init)
    sd_log10_total0 <- rep(0.25, n_grp_init)

    mean_logit_R <- -2
    sd_logit_R <- 1
    sd_bs_R <- 1

    mean_logit_cl <- logit((48 - min_cl) / (max_cl - min_cl))
    sd_logit_cl <- 1
    sd_bs_cl <- 0.5

})


if (names(change) != "none") {
    stan_env[[names(change)]] <- change[[1]]
}



d <- as.list(stan_env)
f <- sampling(plasmofit:::stanmodels$archer_fit, data = as.list(stan_env),
              chains = .n_chains, iter = 1000L + w, warmup = w,
              seed = 538065874,
              save_warmup = FALSE)
write_rds(f, sprintf("wock-fit-%s.rds", names(change)))
cat("Finished fit with", names(change), "altered\n", sep = " ")
cat("\nget_elapsed_time(f):\n")
print(rstan::get_elapsed_time(f))
cat("\n-------------------------------------\n")
cat(names(change), "| leapfrog", round(mean(rstan::get_num_leapfrog_per_iteration(f)), 1),
    "| div", sum(rstan::get_divergent_iterations(f)), "\n")





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

## seconds per gradient evaluation (needs a live DSO)
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

## do divergences concentrate where a scale parameter is small?
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

    list(
        label        = label,
        center_R     = data$center_R,
        center_cl    = data$center_cl,
        fit          = fit,
        lp_by_chain  = sapply(seq_len(ncol(post[, , "lp__"])),
                              function(i) mean(post[, i, "lp__"])),
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

# np <- bayesplot::nuts_params(f)
# div <- np$Iteration[np$Parameter == "divergent__" & np$Value == 1]

# Save fit and dso-dependent output:
res <- summarize_fit(f, d, label = names(change))
write_rds(res, sprintf("wock-fit-RES-%s.rds", names(change)))

