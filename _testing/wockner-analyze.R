suppressPackageStartupMessages({
    library(tidyverse)
    library(rstan)
    library(plasmofit)
    library(bayesplot)  # mcmc_pairs
    library(ggtext)
    library(patchwork)
})


par_names <- c("b_shape", "b_offset", "R", "log10_total0", "cycle_length", "sd_iRBC")


fits <- list.files("_testing", "wock-fit-.*.rds", full.names = TRUE)[c(2,3,1)] |>
    map(read_rds)

posts <- map(fits, \(x) rstan::extract(x, permuted = FALSE))
lps <- map(posts, \(x) sapply(1:4, \(i) mean(x[, i, "lp__"])))


{
    print(do.call(rbind, map(lps, \(x) round(x, 1))))
    print(do.call(rbind, map(posts, \(x) round(apply(x[, , "b_offset[1]"], 2, mean), 3))))
    print(map_dbl(fits, \(x) mean(rstan::get_num_leapfrog_per_iteration(x))))
}


f <- fits[[2]]   # your existing Wockner fit
f0 <- sampling(stan_mod, data = wockner_data, chains = 1, iter = 2,
               warmup = 1, refresh = 0)
dr <- rstan::extract(f, pars = c("b_shape","b_off_vec","log10_total0",
                                 "mu_logit_R","sigma_logit_R","eta_R",
                                 "mu_logit_cl","sigma_logit_cl","eta_cl",
                                 "z_sd_iRBC"), permuted = TRUE)
one <- lapply(dr, \(x) {
    d <- dim(x)
    if (length(d) == 1L) x[1]
    else if (length(d) == 2L) x[1, ]
    else array(x[1, , ], d[-1])
})
u <- rstan::unconstrain_pars(f, one)
system.time(for (i in 1:200) rstan::grad_log_prob(f, u))[["elapsed"]] / 200



{
    rstan::check_treedepth(fits[[2]]) |> print()
    rstan::check_divergences(fits[[2]]) |> print()
    sapply(c("sigma_logit_R","sigma_logit_cl"), \(p) {
        s <- as.matrix(fits[[2]], pars = p)
        c(median = median(s), q05 = quantile(s, .05), q95 = quantile(s, .95))
    }) |> print()
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

map_dbl(fits, fit_cond, max_shape = 250)




