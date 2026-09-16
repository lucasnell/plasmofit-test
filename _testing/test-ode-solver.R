
library(deSolve)
library(rstan)
library(tidyverse)
library(ggtext)


# function parameterized from Kriek et al. 2003
# for sequestration as a function of hpi
gfx = function(age, p1 = 11.3869/467.6209,
               p2 = 1,
               p3 = 18.5802,
               p4 = 0.2242){
    gVal = p1+((p2-p1)/(1+10^(p4*(p3-age))))
    return(gVal)
} # end gfx function

constantPMR.gammaN=function(t, y, parms){
    # Archer model with constant PMR, variable number of age compartments
    cycle_length = parms[1] # in hours
    mu = parms[2] # mortality of non-sequestered iRBCs
    museq = parms[3] # mortality of sequestered iRBCs
    R = parms[4] # parasite multiplication rate (constant)
    n = parms[5] # number of compartments

    lambdaN = n/cycle_length
    lambdaS = lambdaN

    # generate vector of ages through IED:
    ageValues = seq(1/lambdaN,cycle_length,by=1/lambdaN)
    gValues = gfx(ageValues)
    yValues = 1-gValues
    qValues = c(0, 1 - yValues[2:length(yValues)]/yValues[1:(length(yValues)-1)])

    N = y[1:n]
    S = y[(n+1):(2*n)]

    dN = rep(0, n)
    dS = rep(0, n)
    dN[1] = R*(lambdaN * N[n] + lambdaS * S[n]) - (lambdaN + mu) * N[1]
    # Correction from 5/30/24 meeting was
    # dS[1] = - (lambdaS + mu) * S[1]
    # but I believe it should be:
    dS[1] = - (lambdaS + museq) * S[1]
    dN[2:n] = (1-qValues[1:(n-1)]) * lambdaN * N[1:(n-1)] - (lambdaN + mu) * N[2:n]
    dS[2:n] = qValues[1:(n-1)] * lambdaN * N[1:(n-1)] + lambdaS * S[1:(n-1)] - (lambdaS + museq) * S[2:n]

    res=c(dN,dS)
    list(res)
} # end simplified gamma chain model




stan_mod <- stan_model("test-ode.stan")

expose_stan_functions("test-ode.stan")



stan_data = list(nt = 20*24, ts = 1:(20*24),
                 cycle_length = 48,
                 mu = 0,
                 museq = 0,
                 R = 5,
                 n_c = 96L, y0 = rep(0.05, 96L*2L))
parms <- c(stan_data$cycle_length, stan_data$mu, stan_data$museq, stan_data$R,
           stan_data$n_c)




sim_r <- lsoda(stan_data$y0, times = stan_data$ts,
               func = constantPMR.gammaN,
               parms = parms)

sim_st <- lsoda(stan_data$y0, times = stan_data$ts,
                func = \(t, y, parms) {
                    list(constantPMR_gammaN(t, y, parms$cycle_length, parms$mu,
                                            parms$museq, parms$R, parms$n_c))
                },
                parms = stan_data)

all.equal(sim_r, sim_st); max(sim_r - sim_st)



# build_A_R <- function(stan_data) {
#
#     n_c <- stan_data$n_c
#     cycle_length <- stan_data$cycle_length
#     mu <- stan_data$mu
#     museq <- stan_data$museq
#     R <- stan_data$R
#     lambda <- n_c / cycle_length
#
#     q_vals <- (1 - gfx(seq(1/lambda, cycle_length, by=1/lambda))) |>
#         (\(y) c(0, 1 - y[2:length(y)]/y[1:(length(y)-1)]))()
#
#     A <- matrix(0, 2 * n_c, 2 * n_c)
#     for (i in 1:n_c) {
#         A[i, i] = -(lambda + mu)
#         A[n_c+i, n_c+i] = -(lambda + museq)
#     }
#     A[1, n_c] = R * lambda
#     A[1, 2 * n_c] = R * lambda
#     for (i in 2:n_c) {
#         A[i, i-1] = (1 - q_vals[i-1])*lambda
#         A[n_c+i, i-1] = q_vals[i-1]*lambda
#         A[n_c+i, n_c+i-1] = lambda
#     }
#
#     return(A)
# }


# sim_exp_mat <- function(stan_data) {
#
#     if (length(unique(diff(stan_data$ts))) > 1) stop("`sim_exp_mat` expects regular time series")
#
#     q_vals <- with(stan_data,{
#         (1 - gfx(seq(cycle_length/n_c, cycle_length, by = cycle_length/n_c))) |>
#             (\(y) c(0, 1 - y[2:length(y)]/y[1:(length(y)-1)]))()})
#
#     A <- build_A_R(stan_data)
#
#     P <- as.matrix(Matrix::expm(A * 1)) # where dt = 1
#     y <- matrix(0, stan_data$nt, length(stan_data$y0))
#     y[1,] <- stan_data$y0
#     for (t in 2:stan_data$nt) {
#         y[t,] <- P %*% cbind(y[t-1,])
#     }
#
#     return(y)
# }
#
#
# sim_exp_mat_irreg <- function(stan_data) {
#     y0 <- stan_data$y0
#     ts <- stan_data$ts
#     A <- build_A_R(stan_data)
#     y <- matrix(0, length(ts), length(y0))
#     y[1, ] <- y0
#     B <- cbind(y0)                       # 2n x 1, matching Stan's to_matrix()
#     for (i in seq_along(ts)[-1]) {
#         y[i, ] <- as.matrix(Matrix::expm(A * (ts[i] - ts[i-1]))) %*% y[i-1,]
#     }
#     return(y)
# }
#
#
#
# y <- sim_exp_mat(stan_data)
# y2 <- sim_exp_mat_irreg(stan_data)
#
# all.equal(sim_r[,-1], y, check.attributes = FALSE)
# max(abs(sim_r[,-1] - y))
# all.equal(sim_r[,-1], y2, check.attributes = FALSE)
# max(abs(sim_r[,-1] - y2))
#
# all.equal(y, y2, check.attributes = FALSE)
# max(abs(y - y2))



fit <- sampling(stan_mod, data = stan_data,
                algorithm = "Fixed_param", iter = 1, chains = 1)
y_hat <- rstan::extract(fit)$y_hat[1,,]   # nt x 2n

all.equal(sim_r[,-1], y_hat[,], check.attributes = FALSE)



tibble(time = sim_r[,1], circ = rowSums(sim_r[ , 2:(stan_data$n_c + 1)]),
       type = "r") |>
    bind_rows(tibble(time = stan_data$ts, circ = rowSums(y_hat[ , 1:(stan_data$n_c)]),
                     type = "stan")) |>
    ggplot(aes(time, circ, color = type, linetype = type)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(r = "gray60", stan = "dodgerblue"))

tibble(time = sim_r[,1], sequ = rowSums(sim_r[ , (stan_data$n_c + 2):(2*stan_data$n_c + 1)]),
       type = "r") |>
    bind_rows(tibble(time = stan_data$ts,
                     sequ = rowSums(y_hat[ , (stan_data$n_c+1):(2*stan_data$n_c)]),
                     type = "stan")) |>
    ggplot(aes(time, sequ, color = type, linetype = type)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(r = "gray60", stan = "dodgerblue"))



# stan_data2 <- stan_data
# stan_data2[["ts"]] <- sort(sample.int(20*24, 100))
# stan_data2[["nt"]] <- length(stan_data2[["ts"]])
#
# sim_r2 <- lsoda(stan_data2$y0, times = stan_data2$ts,
#                 func = constantPMR.gammaN,
#                 parms = parms)
#
# y2 <- sim_exp_mat_irreg(stan_data2)
#
# all.equal(sim_r2[,-1], y2, check.attributes = FALSE)
# max(abs(sim_r2[,-1] - y2))
#
#
# tibble(time = sim_r2[,1], circ = rowSums(sim_r2[ , 2:(stan_data2$n_c + 1)]),
#        type = "r") |>
#     bind_rows(tibble(time = stan_data2$ts, circ = rowSums(y2[ , 1:(stan_data2$n_c)]),
#                      type = "stan")) |>
#     ggplot(aes(time, circ, color = type)) +
#     geom_line() +
#     scale_color_manual(values = c(r = "gray60", stan = "dodgerblue"))
#
# tibble(time = sim_r2[,1], sequ = rowSums(sim_r2[ , (stan_data2$n_c + 2):(2*stan_data2$n_c + 1)]),
#        type = "r") |>
#     bind_rows(tibble(time = stan_data2$ts,
#                      sequ = rowSums(y2[ , (stan_data2$n_c+1):(2*stan_data2$n_c)]),
#                      type = "stan")) |>
#     ggplot(aes(time, sequ, color = type)) +
#     geom_line() +
#     scale_color_manual(values = c(r = "gray60", stan = "dodgerblue"))




# Test beta_starts ----

betaOffsetStart = function(compartments, shape, offset, initialI){
    times = seq(0, 1, length.out = round(compartments + 1)) + offset
    correctedTimes = times
    correctedTimes[times > 1] = times[times > 1] - 1
    pbetaValues = pbeta(correctedTimes, shape1 = shape, shape2 = shape)
    correctedpbetaValues = pbetaValues
    correctedpbetaValues[times > 1] = pbetaValues[times > 1] + 1
    iVector = initialI * diff(correctedpbetaValues)
    if(round(sum(iVector)) != round(initialI)){
        warning(paste("betaOffsetStart magnitude error (",
                      round(sum(iVector)) - round(initialI), ", offset = ",
                      round(offset, digits = 2), ", shape = ",
                      round(shape, digits = 2), sep = ''), immediate. = T)}
    if(round(length(iVector)) != round(compartments)){
        warning(paste("betaOffsetStart length error, offset = ",
                      round(offset, digits = 2), ", shape = ",
                      shape, sep=''), immediate. = T)}
    return(iVector)
}


for (i in 1:1000) {
    n_c <- sample.int(1000, 1)
    sh <- runif(1)*9+1
    of <- runif(1)
    i0 <- runif(1)*1000
    eq <- isTRUE(all.equal(betaOffsetStart(n_c, sh, of, i0),
                           beta_starts(n_c, sh, of, i0)))
    if (!eq) stop("Noooooo")
}; rm(i, n_c, sh, of, i0, eq)



# fast version (probably not necessary) works well, too
n_c <- 96L
b_shape <- runif(1)*9+1
b_offset <- runif(1)
total0 <- runif(1)*1000
b <- beta_starts(n_c, b_shape, b_offset, total0)
b2 <- beta_starts_fast(n_c, b_shape, b_offset, total0)
all.equal(b, b2)
plot(b, type = "l"); lines(b2, col = "red", lty = 2)




# ============================================================================*
# ============================================================================*
# Testing poly series ----
# ============================================================================*
# ============================================================================*




stan_data2 <- stan_data
stan_data2[["y0"]] <- generate_starts(cycle_length = 48,
                                      n_c = 96L,
                                      b_shape = 5,
                                      b_offset = 0.5,
                                      log10_total0 = log10(100))
stan_data2$nt <- 216L
stan_data2$ts <- stan_data$ts[1:(stan_data2$nt)]


parms2 <- with(stan_data2, { c(cycle_length, mu, museq, R, n_c) })


sim_r2 <- lsoda(stan_data2$y0, times = stan_data2$ts,
               func = constantPMR.gammaN,
               parms = parms2)[,-1]

poly_y <- ew_poly_series(cycle_length = stan_data2$cycle_length, n_c = stan_data2$n_c,
               b_shape = 5, b_offset = 0.5, log10_total0 = log10(100),
               R = stan_data2$R, mu = stan_data2$mu,
               n_obs = stan_data2$nt, ts = stan_data2$ts,
               M = 10, r_max = 10 * 96^2,
               lgam = lgamma(1:(10 * 96^2)),
               r_less1 = (1:(10 * 96^2)) - 1)


tibble(time = stan_data2$ts-1, circ = rowSums(sim_r2[ , 1:(stan_data$n_c)]),
       type = "r") |>
    bind_rows(tibble(time = stan_data2$ts, circ = poly_y, type = "stan")) |>
    ggplot(aes(time, log10(circ), color = type, linetype = type)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(r = "gray60", stan = "dodgerblue")) +
    labs(x = "Time (hours)", y = "log<sub>10</sub>(circulating iRBC)") +
    theme_classic() +
    theme(axis.title.y = element_markdown())



