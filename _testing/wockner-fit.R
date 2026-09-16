# First on local machine:
#
# cd ~/GitHub/Cornell/plasmofit/_testing
# scp wockner-fit.R wockner-cleaned.csv lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/
#
# Now on cluster:
#
# cd /home2/lan68/plasmofit
#
# ... using an interactive job:
# srun -N 1 -n 1 -c 12 --mem=40G --time=24:00:00 --job-name="wock-fit" --pty R --vanilla
#
# ... using a non-interactive job:
#
# cat << EOF > wockner-fit.sh
# #!/bin/bash -l
#
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=12
# #SBATCH --mem=40G
# #SBATCH --time=24:00:00
# #SBATCH --job-name=wock-fit
# #SBATCH --output=wockner-fit.out
# #SBATCH --error=wockner-fit.err
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
# cd ~/GitHub/Cornell/plasmofit/_testing
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/wock-fit-*.rds ./
#



.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_threads)) stop("n_threads cannot be NA")


suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(plasmofit)
})

rstan_options(threads_per_chain = max(1L, n_threads %/% 4L))
options(mc.cores = max(1L, n_threads %/% rstan_options("threads_per_chain")))

stopifnot((options()[["mc.cores"]] * rstan_options("threads_per_chain")) <= n_threads)



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
    center_R <- 1L
    center_cl <- 1L

    # Indices for grouping time series with the same for each:
    grp_init <- summ_by_ts_paras_df$inoc
    grp_R <- summ_by_ts_paras_df$trial
    grp_cl <- summ_by_ts_paras_df$trial
    grp_sd <- summ_by_ts_paras_df$obs_error

    # Parameter limitations:
    max_shape <- 250
    max_R <- 50
    min_cl <- 35
    max_cl <- 55

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
    sd_bs_cl <- 0.1

})



stan_mod <- plasmofit:::stanmodels$archer_fit

# NOTE: THESE TAKE LONG ENOUGH (~10 HRS) THAT THEY SHOULD BE RUN IN SEPARATE JOBS


# for (w in c(400L, 700L, 1000L)) {
#     fit <- sampling(stan_mod, data = as.list(stan_env), chains = 4, refresh = 0,
#                     iter = 1000L + w, warmup = w, save_warmup = FALSE)
#     write_rds(fit, sprintf("wock-fit-%04i.rds", w))
#     cat("Finished fit with", w, "warmup iterations\n", sep = " ")
#     cat("\nget_elapsed_time(fit):\n")
#     print(rstan::get_elapsed_time(fit))
#     cat("\n-------------------------------------")
#     cat("\n-------------------------------------\n")
# }


w = 1000L

fit <- sampling(stan_mod, data = as.list(stan_env), chains = 4,
                iter = 1000L + w, warmup = w, save_warmup = FALSE)
write_rds(fit, sprintf("wock-fit-%04i.rds", w))
cat("Finished fit with", w, "warmup iterations\n", sep = " ")
cat("\nget_elapsed_time(fit):\n")
print(rstan::get_elapsed_time(fit))
cat("\n-------------------------------------")
cat("\n-------------------------------------\n")

