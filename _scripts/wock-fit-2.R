# Combination test: max_cl = 50 AND sd_bs_cl = 0.5 together.
#
# Why: with max_cl = 55, the cycle_length posterior is bimodal (a ~45 h mode and
# a degenerate ~54 h mode pinned at the upper bound). max_cl = 50 alone removed
# that and converged (R-hat 1.005), but left sd_bs_cl = 0.1, a prior tight
# enough that it, not the data, sets the between-trial spread. sd_bs_cl = 0.5
# alone let sigma_logit_cl grow to ~1.14, but only because individual trials
# escaped into the boundary mode (5 of 13 pinned near 55).
#
# This run removes the boundary trap AND gives the between-trial SD room, so
# sigma_logit_cl can finally be read as a statement about the data:
#
#   * sigma_logit_cl stays small, trials stay tightly clustered
#       -> the cycle_length hierarchy is not earning its keep; collapse it to a
#          single shared cycle_length (drops 15 parameters and the funnel).
#   * sigma_logit_cl grows, trials spread out WITHIN (35, 50) without piling up
#     at a bound
#       -> between-trial variation is real; keep the hierarchy.
#
# NOTE, differs deliberately from wockner-fit.R: there, mean_logit_cl was
# computed inside the with() block using max_cl = 55 and the max_cl change was
# applied afterwards, so the max_cl = 50 run actually centred the cycle_length
# prior at 35 + 15 * inv_logit(0.619) = 44.75 h rather than the intended 48 h.
# Here max_cl is set before mean_logit_cl is computed, so the prior is centred
# at 48 h as the expression intends. That makes this a cleaner read on what the
# data want, but it is one more difference from the earlier runs, so keep it in
# mind when comparing.
#
#
# First on local machine:
#
# cd ~/GitHub/Cornell/plasmofit-test
# scp ./_scripts/wock-fit-2.R ./_data/wockner-cleaned.csv lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit/
#
# Now on cluster:
#
# cd /home2/lan68/plasmofit/wock-fit
#
# ... using an interactive job:
# srun -N 1 -n 1 -c 4 --mem=8G --time=24:00:00 --job-name="wock-fit-2" --pty R --vanilla
#
# ... using a non-interactive job (single job, no array needed):
#
# cat << EOF > wock-fit-2.sh
# #!/bin/bash -l
#
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=4
# #SBATCH --mem=8G
# #SBATCH --time=2-00:00:00
# #SBATCH --job-name=wock-fit-2
# #SBATCH --output=wock-fit-2.out
# #SBATCH --error=wock-fit-2.err
# #SBATCH --mail-user=lan68@cornell.edu
# #SBATCH --mail-type=END,FAIL
#
# Rscript --vanilla wock-fit-2.R
#
# EOF
#
# sbatch wock-fit-2.sh
#
# Then when back on local computer:
# cd ~/GitHub/Cornell/plasmofit-test
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit/wock-fit-*.rds ./_data/
#



.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")


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


# Matches the single-change runs, so this is comparable to them:
w <- 700L
.seed <- 538065874

# Both changes, applied together:
.max_cl <- 50
.sd_bs_cl <- 0.5

# Pointwise log likelihood for loo/waic. On here because this config is the
# baseline that simpler model versions (e.g. a single shared cycle_length) get
# compared against, and loo needs log_lik from BOTH sides of the comparison.
# Costs n_total_obs extra numbers per draw, so set to 0L for production fits.
.calc_log_lik <- 1L

.label <- "max_cl50_sd_bs_cl05"


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
    calc_log_lik <- .calc_log_lik
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
    max_cl <- .max_cl        # <- change 1

    # Hyperparameters:
    mean_log_b_shape <- rep(2, n_grp_init)
    sd_log_b_shape <- rep(0.5, n_grp_init)
    mean_log10_total0 <- rep(1, n_grp_init)
    sd_log10_total0 <- rep(0.25, n_grp_init)

    mean_logit_R <- -2
    sd_logit_R <- 1
    sd_bs_R <- 1

    # computed AFTER max_cl is set, so this really is centred at 48 h
    mean_logit_cl <- logit((48 - min_cl) / (max_cl - min_cl))
    sd_logit_cl <- 1
    sd_bs_cl <- .sd_bs_cl    # <- change 2

})


d <- as.list(stan_env)

stopifnot(d$max_cl == 50, d$sd_bs_cl == 0.5)
cat("prior for cycle_length centred at",
    round(d$min_cl + (d$max_cl - d$min_cl) * plogis(d$mean_logit_cl), 2), "h\n")


f <- sampling(plasmofit:::stanmodels$archer_fit, data = d,
              chains = .n_chains, iter = 1000L + w, warmup = w,
              seed = .seed,
              save_warmup = FALSE)
write_rds(f, sprintf("wock-fit-%s.rds", .label))
cat("Finished fit with", .label, "\n", sep = " ")
cat("\nget_elapsed_time(f):\n")
print(rstan::get_elapsed_time(f))
cat("\n-------------------------------------\n")
cat(.label, "| leapfrog", round(mean(rstan::get_num_leapfrog_per_iteration(f)), 1),
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

# Save fit and dso-dependent output:
res <- summarize_fit(f, d, label = .label)
write_rds(res, sprintf("wock-fit-RES-%s.rds", .label))




## ---- the readout this run exists for -------------------------------------

cat("\n=====================================================\n")
cat("max R-hat   :", round(max(res$diagnostics$rhat, na.rm = TRUE), 4), "\n")
cat("divergences :", res$n_divergent, "\n")
cat("cond        :", round(res$cond, 2), "\n")
cat("lp by chain :", paste(round(res$lp_by_chain, 1), collapse = "  "), "\n")

cat("\n-- sigma_logit_cl: posterior vs prior --\n")
s <- as.matrix(f, pars = "sigma_logit_cl")
cat(sprintf("  posterior median %.4f  (5%%-95%%: %.4f - %.4f)\n",
            median(s), quantile(s, .05), quantile(s, .95)))
set.seed(1)
pr <- abs(rnorm(4e4, 0, d$sd_bs_cl))
cat(sprintf("  prior     median %.4f  (5%%-95%%: %.4f - %.4f)\n",
            median(pr), quantile(pr, .05), quantile(pr, .95)))
cat("  -> posterior well below prior means the data do NOT want between-trial\n")
cat("     variation, and the cycle_length hierarchy can be collapsed.\n")

cat("\n-- per-trial cycle_length (watch for piling up at the bounds) --\n")
cl <- as.matrix(f, pars = "cycle_length")
q <- apply(cl, 2, quantile, c(.05, .5, .95))
for (j in seq_len(ncol(cl))) {
    flag <- if (q[2, j] > d$max_cl - 1.5) "   <-- near max_cl" else
            if (q[2, j] < d$min_cl + 1.5) "   <-- near min_cl" else ""
    cat(sprintf("  trial %2d: %5.2f  (%5.2f - %5.2f)%s\n",
                j, q[2, j], q[1, j], q[3, j], flag))
}
cat(sprintf("\n  spread of per-trial medians: %.2f h\n", diff(range(q[2, ]))))

## elpd for comparison against simpler model versions later.
## Save the loo object, not just the number: loo::loo_compare() needs the
## objects, and recomputing it needs log_lik which we may not keep around.
if (.calc_log_lik == 1L && requireNamespace("loo", quietly = TRUE)) {
    ll <- loo::extract_log_lik(f, parameter_name = "log_lik", merge_chains = FALSE)
    reff <- loo::relative_eff(exp(ll))
    lll <- loo::loo(ll, r_eff = reff)
    write_rds(lll, sprintf("wock-fit-LOO-%s.rds", .label))
    cat("\n-- loo --\n")
    print(lll)
    cat("\n  (compare later with: loo::loo_compare(loo_A, loo_B))\n")
    cat("  any Pareto k > 0.7:", sum(loo::pareto_k_values(lll) > 0.7), "of", ncol(ll), "\n")
} else if (.calc_log_lik == 1L) {
    cat("\n  [loo package not available; log_lik is saved in the fit]\n")
}
cat("=====================================================\n")
