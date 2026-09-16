
suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(bayesplot)  # mcmc_pairs
    library(ggtext)
    library(patchwork)
})


options(mc.cores = max(1L, parallel::detectCores()-2L))

stan_mod <- stan_model("archer-fit-single.stan")


# Time points sampled from Wockner:
ts <- seq(72, 216, 12)


# # Just expose functions in archer-functions.stan
# expose_stan_functions("test-ode.stan")


build_A_R <- function(cycle_length, mu, museq, q_vals, R, n_c) {

    lambda <- n_c / cycle_length

    A <- matrix(0, 2 * n_c, 2 * n_c)
    for (i in 1:n_c) {
        A[i, i] = -(lambda + mu)
        A[n_c+i, n_c+i] = -(lambda + museq)
    }
    A[1, n_c] = R * lambda
    A[1, 2 * n_c] = R * lambda
    for (i in 2:n_c) {
        A[i, i-1] = (1 - q_vals[i-1])*lambda
        A[n_c+i, i-1] = q_vals[i-1]*lambda
        A[n_c+i, n_c+i-1] = lambda
    }

    return(A)
}

y_vals_R <- function(cycle_length, n_c) {
    ages <- seq(cycle_length/n_c, cycle_length, length.out = n_c)
    p1 <- 11.3869/467.6209; p2 <- 1; p3 <- 18.5802; p4 <- 0.2242
    1 - (p1 + (p2 - p1)/(1 + 10^(p4*(p3 - ages))))
}

beta_starts_R <- function(n_c, b_shape, b_offset, total0) {
    times <- seq(0, 1, length.out = n_c + 1) + b_offset
    ct <- ifelse(times > 1, times - 1, times)
    v  <- pbeta(ct, b_shape, b_shape) + ifelse(times > 1, 1, 0)
    total0 * diff(v)
}

generate_starts_R <- function(cycle_length, n_c, b_shape, b_offset, total0) {
    yv <- y_vals_R(cycle_length, n_c)
    w  <- beta_starts_R(n_c, b_shape, b_offset, total0)
    c(yv * w, (1 - yv) * w)
}

q_vals_R <- function(cycle_length, n_c) {
    yv <- y_vals_R(cycle_length, n_c)
    c(0, 1 - yv[-1]/yv[-n_c])
}

simulate_circ <- function(truth, n_c, ts, dt) {

    # n_c = .n_c; ts = seq(3,9,0.5) |> round(1); dt = 0.5
    # rm(n_c, ts, dt, q, A, P, n_lik, idx, s, out, k)

    cycle_length <- truth$cycle_length

    q  <- q_vals_R(cycle_length, n_c)
    A  <- build_A_R(cycle_length, 0, 0, q, truth$R, n_c)   # your verified builder
    P  <- as.matrix(Matrix::expm(A * dt))
    n_lik <- round(ts[length(ts)]/dt) + 1
    idx   <- round(ts/dt) + 1
    s <- generate_starts_R(cycle_length, n_c, truth$b_shape, truth$b_offset, truth$total0)
    out <- numeric(n_lik)
    out[1] <- sum(s[1:n_c])
    for (k in 2:n_lik) {
        s <- P %*% s
        out[k] <- sum(s[1:n_c])
    }
    return(out[idx])
}

make_coef <- function(cycle_length, n_c, ts, M) {
    lambda <- n_c / cycle_length
    q  <- q_vals_R(cycle_length, n_c)          # published convention, q[1] = 0
    ck <- c(1, cumprod(1 - q[1:(n_c - 1)]))    # ck[k] = prod_{j<k}(1 - q_j)

    lp <- function(r, lt) (r - 1) * log(lt) - lt - lgamma(r)

    coef <- array(0, dim = c(length(ts), M + 1, 2 * n_c))
    for (i in seq_along(ts)) {
        stopifnot(ts[i] > 0)                   # t = 0 needs special handling
        lt <- lambda * ts[i]
        for (j in 1:n_c) {                     # m = 0: initial circulating cohort
            ks <- j:n_c
            coef[i, 1, j] <- sum((ck[ks] / ck[j]) * exp(lp(ks - j + 1, lt)))
        }
        ks <- 1:n_c
        for (m in 1:M) for (j in 1:n_c) {      # m >= 1: after m bursting events
            v <- sum(ck[ks] * exp(lp(m * n_c - j + ks + 1, lt)))
            coef[i, m + 1, j]       <- v       # from initial N_j
            coef[i, m + 1, n_c + j] <- v       # from initial S_j (same transit)
        }
    }
    coef
}



# To plot prior distributions
# plot_prior_distr <- function(m, s, xmax = 1, xmin = 0) {
#     x <- seq(xmin, xmax, length.out = .n_c)
#     y <- dnorm(x, m, s)
#     y <- y / sum(y)
#     plot(x, y, type = "l", ylim = c(0, max(y)))
# }
#
#
# d <- sapply(1:1000, \(i) {
#     bs <- runif(1) * 9 + 1
#     bo <- runif(1)
#     t0 <- runif(1, 1, 10000)
#     cl <- runif(1,  24, 72)
#     nc <- sample(50:300, 1)
#     x <- generate_starts(cycle_length = cl, n_c = nc, b_shape = bs, b_offset = bo, total0 = t0)
#     y <- generate_starts_R(cycle_length = cl, n_c = nc, b_shape = bs, b_offset = bo, total0 = t0)
#     # if (!isTRUE(all.equal(x, y, tolerance = 1e-4))) stop("Nooooooo")
#     mean(pmin(abs(x-y), abs(x)))
# })
# max(d)




# with n_c = 48:
#   - Took ~35 sec
#   - Weird parameter estimates - do not use
#
# with n_c = 96:
#   - Took ~60 sec
#   - Good parameter estimates
#   - Passes all diagnostics
#
# with n_c = 192:
#   - Took ~250 sec (~640 with adapt_delta = 0.99)
#   - Okay parameter estimates
#   - Fails diagnostics (divergent outcomes, rhat, ess)
#   - 0.75% divergence with adapt_delta = 0.99, others pass
#   - Interestingly, contractions are higher :
#       b_shape     b_offset log10_total0            R cycle_length
#     0.6566489    0.9892901    0.5423014    0.9799984    0.9729279
#
.n_c <- 48L * 2L

truth <- list(b_shape = 100, b_offset = 0.30, total0 = 1e4, R = 6.4,
              cycle_length = 48, sd_iRBC = 0.15) |>
    (\(x) {x$log10_total0 <- log10(x$total0); x} )()
y_true <- simulate_circ(truth, n_c = .n_c, ts = ts, dt = 12)
y_clean <- y_true                      # no noise

# set.seed(42)
y_obs <- pmax(10^(log10(y_true + 1) + rnorm(length(y_true), 0, truth$sd_iRBC)) - 1, 0)


par_names <- c("b_shape", "b_offset", "log10_total0", "R", "cycle_length", "sd_iRBC")


stan_data_sim <- list(
    n_obs = length(y_obs),             # Number of observed time points
    ts = ts,  # Times to evaluate (in hours)
    mu = 0,       # mortality of all iRBCs (set to 0)
    n_c = .n_c,               # Number of compartments
    y = y_obs,       # iRBC abundance (untransformed)
    dt_full = 1,  # Step size for full output
    max_shape = 250,  # set upper limit on beta shape parameter
    t_ref = 120L,
    # ========== Hyperparameters ==========
    mu_shape = 150,  # mean of shape parameter
    sd_shape = 100,  # stdev of shape parameter
    mu_offset = 0.5,  # mean of beta offset parameter
    sd_offset = 0.5,  # stdev of beta offset parameter
    mu_total = 4,  # mean of log10(total starting abundance)
    sd_total = 0.2,  # stdev of log10(total starting abundance)
    mu_R = 5,  # mean of parasite multiplication rate
    sd_R = 5,  # stdev of parasite multiplication rate
    mu_cl = 45,  # mean of cycle length
    sd_cl = 5,  # stdev of cycle length
    min_cl = 35,  # minimum cycle length (required for `M` below)
    max_cl = 55  # maximum cycle length (required to prevent `y` overflow)
)


# Takes ~ 1 min
fit <- sampling(stan_mod, data = stan_data_sim, chains = 4, iter = 2000,
                # seed = 1854403876,
                init = function() list(b_shape = runif(1, 20, 250),
                                       b_offset = runif(1, 0, 1),
                                       log10_total0 = runif(1, 3, 5),
                                       R = runif(1, 3, 10)))
# fit


# cor(as.numeric(extract(fit, "cycle_length", permute = FALSE)[,,1]), as.numeric(extract(fit, "phi_ref", permute = FALSE)[,,1]))
# # [1] -0.9305012


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
        mcmc_areas(fit, p) + geom_vline(xintercept = truth[[p]], linetype = "22")
    }) |>
    wrap_plots(ncol = 2)



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
tibble(par = par_names,
       diags = map(par, \(p) {
           sims <- extract(fit, p, permuted = FALSE)[,,1]
           rh <- Rhat(sims)
           ess_b <- ess_bulk(sims)
           ess_t <- ess_tail(sims)
           tibble(rhat = rh, ess_bulk = ess_b, ess_tail = ess_t)
       })) |>
    unnest(diags) |>
    mutate(rhat = num(rhat, digits = 4))
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

