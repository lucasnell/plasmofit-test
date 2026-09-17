
# First on local machine:
#
# cd ~/GitHub/Cornell/plasmofit/_testing
# scp cluster-test-archer-fit.R lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/plasmofit/
#
# Now on cluster:
#
# cd /home2/lan68/plasmofit
#
# ... using an interactive job:
# srun -N 1 -n 1 -c 12 --mem=40G --time=24:00:00 --job-name="plasmofit" --pty R --vanilla
#
# ... using a non-interactive job:
#
# cat << EOF > plasmofit.sh
# #!/bin/bash -l
#
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=12
# #SBATCH --mem=40G
# #SBATCH --time=24:00:00
# #SBATCH --job-name=plasmofit
# #SBATCH --output=plasmofit.out
# #SBATCH --error=plasmofit.err
# #SBATCH --mail-user=lan68@cornell.edu
# #SBATCH --mail-type=END,FAIL
#
# Rscript --vanilla cluster-test-archer-fit.R
#
# EOF
#
# sbatch plasmofit.sh
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



# Time points sampled from Wockner:
all_ts <- seq(72, 216, 12)


simulate_circ <- function(n_c, cycle_length, b_shape, b_offset, log10_total0, R,
                          ts, dt) {

    # If you ever want to do this in R, use below. I want to pass it in this
    # function because it operates on only a single time series.
    # gcd <- function(x, y) {
    #     ifelse(y == 0, x, gcd(y, x %% y))
    # }
    # iters <- 0L
    # while (!all.equal(round(ts), ts)) {
    #     ts <- ts * 1e3
    #     iters <- iters + 1L
    #     if (iters > 3L) stop("times cannot be converted to integers")
    # }
    # dt <- Reduce(gcd, ts)

    y0 <- plasmofit:::generate_starts(cycle_length, n_c, b_shape, b_offset,
                          log10_total0)
    plasmofit:::mat_exp_series(y0 = y0, cycle_length = cycle_length, n_c = n_c,
                   R = R, mu = 0, ts = ts, dt = dt)
}





# =============================================================================*
# =============================================================================*
# Setup simulations
# =============================================================================*
# =============================================================================*


stan_env <- new.env()

with(stan_env, {

    set.seed(2104249041)

    # Numbers of groups:
    n_ts <- 64L
    n_grp_init <- 2L
    n_grp_R <- 8L
    n_grp_cl <- 8L
    n_grp_sd <- 1L

    # Other integers:
    n_obs <- rep(8L, n_ts)
    ts <- map(n_obs, \(n) sort(sample(all_ts, n))) |> list_c()
    n_total_obs <- sum(n_obs)
    mu <- 0
    n_c <- 96L
    dt_full <- 12
    run_check <- 0L
    center_R <- 1L
    center_cl <- 1L

    # Indices for grouping time series with the same for each:
    grp_init <- rep(1:n_grp_init, each = n_ts / n_grp_init)
    grp_R <- rep(1:n_grp_R, n_ts / n_grp_R)
    grp_cl <- rep(rep(1:n_grp_cl, each = 2L), n_ts / (2*n_grp_cl))
    grp_sd <- rep(1, n_ts)

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

# Actual true parameters

# not sure why, but I was getting weird errors if creating copy of stan_env
# directly
truth_env <- new.env()
for (name in ls(stan_env, all.names = TRUE)) {
    assign(name, get(name, stan_env), truth_env)
}; rm(name)

with(truth_env, {

    set.seed(1333619777)

    b_shape <- rlnorm(n_grp_init, mean_log_b_shape, sd_log_b_shape)
    b_offset <- runif(n_grp_init, 0.1, 0.9)
    log10_total0 <- rnorm(n_grp_init, mean_log10_total0, sd_log10_total0)

    mu_logit_R <- rnorm(1, mean_logit_R, sd_logit_R)
    sigma_logit_R <- rtnorm(1, lower = 0, sigma = sd_bs_R)
    z_R  <- rnorm(n_grp_R)

    mu_logit_cl  <- rnorm(1, mean_logit_cl, sd_logit_cl)
    sigma_logit_cl <- rtnorm(1, lower = 0, sigma = sd_bs_cl)
    z_cl  <- rnorm(n_grp_cl)

    R  <- max_R * inv_logit(mu_logit_R + z_R * sigma_logit_R)
    cycle_length <- min_cl + (max_cl - min_cl) * inv_logit(mu_logit_cl + z_cl * sigma_logit_cl)

    sd_iRBC <- rtnorm(n_grp_sd, lower = 0, sigma = 0.2)

    start_idx <- c(1, cumsum(head(n_obs, -1))+1)

    y_true <- numeric(n_total_obs)

    for (i in 1:n_ts) {
        s <- start_idx[i]
        e <- s + n_obs[i] - 1L
        bs_i <- b_shape[grp_init[i]]
        bo_i <- b_offset[grp_init[i]]
        lt_i <- log10_total0[grp_init[i]]
        cl_i <- cycle_length[grp_cl[i]]
        R_i <- R[grp_R[i]]
        y_true[s:e] <- simulate_circ(n_c, cl_i, bs_i, bo_i, lt_i, R_i,
                                     ts[s:e], dt_full)
    }; rm(i, s, e, bs_i, bo_i, lt_i, cl_i, R_i)

})

set.seed(831910407)
stan_env$y <- pmax(0, 10^(log10(truth_env$y_true + 1) +
                           rnorm(stan_env$n_total_obs, 0, truth_env$sd_iRBC)) - 1)


stan_mod <- plasmofit:::stanmodels$archer_fit



fit6 <- sampling(stan_mod,
                 data = as.list(stan_env),
                 chains = 4, iter = 1300, cores = 4, warmup = 300, save_warmup = FALSE)
write_rds(fit6, "fit6.rds")

fit7 <- sampling(stan_mod,
                 data = as.list(stan_env),
                 chains = 4, iter = 1400, warmup = 1000, save_warmup = FALSE)
write_rds(fit7, "fit7.rds")

fit8 <- sampling(stan_mod,
                 data = as.list(stan_env),
                 chains = 4, iter = 2000, warmup = 1000, save_warmup = FALSE)
write_rds(fit8, "fit8.rds")









fits <- map(6:8, \(x) read_rds(sprintf("_testing/fit%i.rds", x)))

# Alternative to above:
fits <- map2(i = c(1300, 1400, 2000), w = c(300, 1000, 1000),
             \(i, w) {
                 sampling(stan_mod, data = as.list(stan_env),
                          chains = 4, iter = i, warmup = w, save_warmup = FALSE)
             })

posts <- map(fits, \(x) rstan::extract(x, permuted = FALSE))

lps <- map(posts, \(x) sapply(1:4, \(i) mean(x[, i, "lp__"])))

{
    print(map(lps, \(x) round(x, 1)))
    print(map(posts, \(x) round(apply(x[, , "b_offset[1]"], 2, mean), 3)))
    print(map(posts, \(x) round(apply(x[, , "b_offset[2]"], 2, mean), 3)))
}


par_names <- c("b_shape", "b_offset", "R", "log10_total0", "cycle_length", "sd_iRBC")

for (i in 1:length(fits)) {
    par_names |>
        map(\(p) {
            sims <- extract(fits[[i]], p, permuted = FALSE)
            rh <- apply(sims, 3, Rhat)
            ess_b <- apply(sims, 3, ess_bulk)
            ess_t <- apply(sims, 3, ess_tail)
            tibble(par = names(rh), rhat = rh, ess_bulk = ess_b, ess_tail = ess_t)
        }) |>
        list_rbind() |>
        mutate(rhat = num(rhat, digits = 4)) |>
        print(n = 25)
}



fit_cond <- function(f, max_shape) {
    post <- as.matrix(f, pars = c("b_shape","b_offset","log10_total0",
                                  "mu_logit_R","sigma_logit_R","eta_R",
                                  "mu_logit_cl","sigma_logit_cl","eta_cl",
                                  "z_sd_iRBC"))
    lg <- plasmofit::logit
    j <- grep("^b_shape",  colnames(post)); post[, j] <- lg((post[, j] - 2) / (max_shape - 2))
    j <- grep("^b_offset", colnames(post)); post[, j] <- lg(post[, j])
    j <- grep("^sigma_",   colnames(post)); post[, j] <- log(post[, j])
    ev <- eigen(cor(post), only.values = TRUE)$values
    sqrt(max(ev) / min(ev))
}

map(fits, fit_cond, max_shape = stan_env$max_shape)




