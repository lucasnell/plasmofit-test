
suppressPackageStartupMessages({
    library(plasmofit)
    library(tidyverse)
    library(rstan)
    library(bayesplot)  # mcmc_pairs
    library(ggtext)
    library(patchwork)
})


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

    # Indices for grouping time series with the same for each:
    grp_init <- rep(1:n_grp_init, each = n_ts / n_grp_init)
    grp_R <- rep(1:n_grp_R, n_ts / n_grp_R)
    grp_cl <- rep(rep(1:n_grp_cl, each = 2L), n_ts / (2*n_grp_cl))
    grp_sd <- rep(1, n_ts)

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


# timing_run <- function(data, iter = 60, ...) {
#     f <- sampling(stan_mod, data = data, chains = 1,
#                   iter = iter, warmup = iter/2, save_warmup = FALSE,
#                   refresh = 0,
#                   control = list(adapt_engaged = FALSE, stepsize = 0.05), ...)
#     nl <- rstan::get_num_leapfrog_per_iteration(f)
#     c(mean_leapfrog = mean(nl),
#       sec_per_grad  = sum(rstan::get_elapsed_time(f)[,"sample"]) / sum(nl))
# }
#
#
# data = modifyList(as.list(stan_env), list(center_R = 1L, center_cl = 1L))
# iter = 60
#
# mean(rstan::get_sampler_params(f, inc_warmup = FALSE)[[1]][,"stepsize__"])
#
# f <- sampling(stan_mod, data = data, chains = 1,
#               iter = iter, warmup = iter/2, save_warmup = FALSE,
#               refresh = 0,
#               control = list(adapt_engaged = FALSE, stepsize = 0.001))
# rstan::get_divergent_iterations(f) |> mean()
#
# nl <- rstan::get_num_leapfrog_per_iteration(f)
# c(mean_leapfrog = mean(nl),
#   sec_per_grad  = sum(rstan::get_elapsed_time(f)[,"sample"]) / sum(nl))
#
#
#
# timing_run(modifyList(as.list(stan_env), list(center_R = 1L, center_cl = 1L)))
# timing_run(modifyList(as.list(stan_env), list(center_R = 1L, center_cl = 0L)))
# timing_run(modifyList(as.list(stan_env), list(center_R = 0L, center_cl = 1L)))
# timing_run(modifyList(as.list(stan_env), list(center_R = 0L, center_cl = 0L)))
#
#
#
# # # Test with 1 chain, 400 iterations took ~46 min
# # # (10.7 s/iteration in warmup and 14.0 s/iteration in sampling)
# # #
# # #
# # # fit <- sampling(stan_mod, data = modifyList(as.list(stan_env), list(run_check = 0)),
# # #                 chains = 1, iter = 20)
# # fit0 <- sampling(stan_mod, data = as.list(stan_env),
# #                  chains = 1, iter = 20, save_warmup = FALSE)
# # rstan::get_elapsed_time(fit0)



# # Took ~38 min
# fit <- sampling(stan_mod,
#                  data = modifyList(as.list(stan_env),
#                                    list(center_R = 1L, center_cl = 1L)),
#                  chains = 4, iter = 1300, cores = 4, warmup = 300, save_warmup = FALSE)
#
# write_rds(fit, "_testing/fit4.rds")

fit <- read_rds("_testing/fit4.rds")
fit4 = fit

# fit



par_names <- c("b_shape", "b_offset", "R", "log10_total0", "cycle_length", "sd_iRBC")


# # cor(as.numeric(extract(fit1, "cycle_length", permute = FALSE)[,,1]), as.numeric(extract(fit1, "b_offset", permute = FALSE)[,,1]))
# # # [1] -0.5046214
#
# post <- as.matrix(fit, pars = c("b_shape","b_offset","log10_total0",
#                                  "mu_logit_R","sigma_logit_R","z_R",
#                                  "mu_logit_cl","sigma_logit_cl","z_cl","z_sd_iRBC"))
# post[,"b_offset[1]"] <- qlogis(post[,"b_offset[1]"])   # and [2]
# post[,"b_shape[1]"]  <- qlogis((post[,"b_shape[1]"] - 2) / (stan_env$max_shape - 2))
# post[,grep("sigma_", colnames(post))] <- log(post[,grep("sigma_", colnames(post))])
# ev <- eigen(cor(post), only.values = TRUE)$values
# sqrt(max(ev) / min(ev))
#
# e <- eigen(cor(post))
# round(e$vectors[, 27][order(abs(e$vectors[, 27]), decreasing = TRUE)][1:5], 3)
# colnames(post)[order(abs(e$vectors[, 27]), decreasing = TRUE)[1:5]]














# The method works well to approximate the exponential approach to the ODE model
# Anything above ~ 1e-6 is a potential problem
extract(fit, "max_rel_diff")[[1]] |> max()
# [1] 5.173462e-13


# This should be added to diagnostics to see if offset should be rotated
q <- quantile(as.matrix(fit, pars = "b_offset"), c(0.005, 0.995))
if (q[1] < 0.02 || q[2] > 0.98) warning("b_offset posterior crowds the boundary; rotate the interval")


# extract(fit, "cycle_length")[[1]] |> hist()

# bayesplot::mcmc_trace(fit, "cycle_length")
map(par_names, \(p) {
        mcmc_areas(fit, p) + geom_vline(xintercept = truth_env[[p]], linetype = "22")
    }) |>
    wrap_plots(ncol = 2)


{
    post <- rstan::extract(fit, permuted = FALSE)
    print(round(apply(post[, , "b_offset[1]"], 2, mean), 3))
    print(round(apply(post[, , "cycle_length[1]"], 2, mean), 3))
    print(round(apply(post[, , "R[1]"], 2, mean), 3))
}



# Test for multimodality:
tibble(chain = 1:4) |>
    (\(d) {
        for (p in par_names) {
            d[[p]] <- extract(fit, p, permuted = FALSE)[,,1] |>
                colMeans() |> as.numeric()
        }
        return(d)
    })()

# Looks like no multimodality!


check_divergences(fit); cat("\n") # want 0
# 0 of 4000 iterations ended with a divergence.

get_bfmi(fit)  # want > 0.3
# [1] 0.8146974 0.8302787 0.7932742 0.9061616

#
# Ideally Rhat < 1.01, bulk and tail ESS > 400
#
par_names |>
    map(\(p) {
           sims <- extract(fit, p, permuted = FALSE)
           rh <- apply(sims, 3, Rhat)
           ess_b <- apply(sims, 3, ess_bulk)
           ess_t <- apply(sims, 3, ess_tail)
           tibble(par = names(rh), rhat = rh, ess_bulk = ess_b, ess_tail = ess_t)
       }) |>
    list_rbind() |>
    mutate(rhat = num(rhat, digits = 4)) |>
    print(n = 25)
# # A tibble: 6 × 4
#   par               rhat ess_bulk ess_tail
#   <chr>        <num:.4!>    <dbl>    <dbl>
# 1 b_shape         1.0006    2201.    1538.
# 2 b_offset        1.0032     865.     780.
# 3 R               1.0021    1132.    1017.
# 4 log10_total0    1.0021    1443.    1733.
# 5 cycle_length    1.0030     939.     880.
# 6 sd_iRBC         1.0009    1804.    1471.


post <- as.data.frame(fit, pars = par_names)

# posterior z-score: how far the posterior sits from the truth
z <- sapply(par_names, function(p)
    (mean(post[[p]]) - truth[[p]]) / sd(post[[p]]))

abs(z)

# > ... the posterior z-score for each parameter... measures how closely the
# > posterior recovers the parameters of the true data generating process...
# > The smaller the absolute z-score, the closer the bulk of the posterior is
# > to the true parameter: z-scores beyond the absolute value of three to four
# > may indicate substantial bias (Schad et al., 2021).
# > ...
# > We found that most of our model parameters... and hyperpriors... were
# > estimated with accuracy, precision and identifiability, with the absolute
# > posterior z-scores well below three... The strain-level variation... tended
# > towards overfitting with the posterior z-score of 2:99...
# > Thus, caution might be warranted when interpreting strain-specific
# > differences in this parameter.
#
# Kamiya et al. (2021)



# posterior contraction: how much the data added over the prior
prior_var <- with(stan_data_sim, {c(b_shape = sd_shape^2,
                                    b_offset = sd_offset^2,
                                    R = sd_R^2,
                                    log10_total0 = sd_total^2,
                                    cycle_length = sd_cl^2,
                                    sd_iRBC = NA)})
contraction <- 1 - sapply(par_names, function(p) var(post[[p]]) / prior_var[[p]])
contraction

# > The posterior contraction values close to zero indicate that data contain
# > little information (i.e., poor identifiability, rendering priors strongly
# > informative). Conversely, values close to one indicate that data are much
# > more informative than the prior (Schad et al., 2021).
# > ...
# > We found that most of our model parameters... and hyperpriors... were
# > estimated with accuracy, precision and identifiability, with ... posterior
# > contraction values beyond 75%... One parameter... showed a comparatively
# > lower posterior contraction value, yet its posterior distribution
# > contracted by 64.0%, meaning that data still provided substantial
# > information over the prior distribution
# Kamiya et al. (2021)


mcmc_pairs(fit, c("log10_total0", "R"))

mcmc_pairs(fit, c("b_offset", "cycle_length"))


mcmc_acf_bar(fit, "R", lags = 5)





fit_df <- tibble(time = seq(0, max(ts), stan_data_sim$dt_full),
       y_pred = extract(fit, "y_hat_full")[[1]] |> apply(2, median),
       y_lo = extract(fit, "y_hat_full")[[1]] |> apply(2, quantile, probs = 0.025),
       y_hi = extract(fit, "y_hat_full")[[1]] |> apply(2, quantile, probs = 0.975))

tibble(time = ts, y = log10(1 + y_obs)) |>
    ggplot(aes(time, y)) +
    geom_point(shape = 19, size = 1.5) +
    geom_line(data = fit_df, aes(x=time, y=log10(1+y_pred)), color = "red") +
    geom_line(data = fit_df, aes(x=time, y=log10(1+y_lo)), color = "red", linetype=3) +
    geom_line(data = fit_df, aes(x=time, y=log10(1+y_hi)), color = "red", linetype=3) +
    labs(x = "Time (hrs)", y = "log<sub>10</sub>(1 + iRBCs)") +
    theme_classic() +
    theme(axis.title.y = element_markdown())




# ===================================================================*
# ===================================================================*
# Fit real dataset ----
# ===================================================================*
# ===================================================================*

source("_testing/Wockner-read.R")  # creates paras_df


best_df <- paras_df |>
    group_by(id) |>
    mutate(n_obs = sum(!is.na(para)),
           first_obs = time[which(!is.na(para))[1]]) |>
    ungroup() |>
    filter(n_obs == max(n_obs)) |>
    filter(first_obs == min(first_obs)) |>
    filter(id == id[1]) |>
    filter(!is.na(para))

best_ts <- best_df[["time"]]
best_y <- best_df[["para"]]

best_data <- stan_data_sim
best_data$mu_total <- 0
best_data$sd_total <- 1
best_data$n_obs = length(best_y)
best_data$ts = best_ts
best_data$y <- best_y
best_data$t_ref <- best_ts[which(abs(best_ts - 120L) == min(abs(best_ts - 120L)))[1]]
# Takes ~ 1 min
best_fit <- sampling(stan_mod, data = best_data, chains = 4, iter = 2000,
                seed = 1854403876,
                init = function() list(b_shape = runif(1, 20, 250),
                                       b_offset = runif(1, 0, 1),
                                       log10_total0 = runif(1, 3, 5),
                                       R = runif(1, 3, 10)))
# best_fit


extract(best_fit, "max_rel_diff")[[1]] |> max()


map(par_names, \(p) mcmc_areas(best_fit, p)) |>
    wrap_plots(ncol = 2)



# Test for multimodality:
tibble(chain = 1:4) |>
    (\(d) {
        for (p in par_names) {
            d[[p]] <- extract(best_fit, p, permuted = FALSE)[,,1] |>
                colMeans() |> as.numeric()
        }
        return(d)
    })()



check_divergences(best_fit); cat("\n") # want 0

get_bfmi(best_fit)  # want > 0.3


#
# Ideally Rhat < 1.01, bulk and tail ESS > 400
#
tibble(par = par_names,
       diags = map(par, \(p) {
           sims <- extract(best_fit, p, permuted = FALSE)[,,1]
           rh <- Rhat(sims)
           ess_b <- ess_bulk(sims)
           ess_t <- ess_tail(sims)
           tibble(rhat = rh, ess_bulk = ess_b, ess_tail = ess_t)
       })) |>
    unnest(diags) |>
    mutate(rhat = num(rhat, digits = 4))



best_post <- as.data.frame(best_fit, pars = par_names)

# # posterior z-score: how far the posterior sits from the truth
# z <- sapply(par_names, function(p)
#     (mean(post[[p]]) - truth[[p]]) / sd(post[[p]]))
#
# abs(z)


# posterior contraction: how much the data added over the prior
best_prior_var <- with(best_data, {c(b_shape = sd_shape^2,
                                    b_offset = sd_offset^2,
                                    R = sd_R^2,
                                    log10_total0 = sd_total^2,
                                    cycle_length = sd_cl^2,
                                    sd_iRBC = NA)})
best_contraction <- 1 - sapply(par_names, function(p) var(post[[p]]) / best_prior_var[[p]])
best_contraction

