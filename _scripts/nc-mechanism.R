## WHY does cycle_length shift with n_c? A deterministic test.
##
## Everything here is noiseless and has no priors, no sampler and no
## estimation: it asks what the FORWARD MAP does. Generate a trajectory at one
## n_c with a known cycle_length, then find the cycle_length that best
## reproduces it at another n_c. Whatever shift appears is built into the
## model's arithmetic, not into inference.
##
## Two candidate channels, both of which n_c controls:
##
##  (1) DESYNCHRONISATION. Transit over a cycle is Erlang(n_c, n_c/cl), so the
##      stage distribution's sd after k cycles is cl*sqrt(k/n_c). Raising n_c
##      slows the decay of the oscillation. This is a genuine biological
##      assumption and does NOT vanish as n_c grows -- at n_c = Inf the chain
##      is deterministic and there is no desynchronisation at all.
##
##  (2) SEQUESTRATION-WINDOW DISCRETISATION. make_log_y_vals puts the
##      circulating fraction on a logistic in ABSOLUTE developmental age,
##      centred at p3 = 18.5802 h, and evaluates it at age k*cl/n_c for
##      compartment k. The grid spacing is cl/n_c hours, and the Archer
##      convention q[1] = 0 adds a documented O(1/n_c) delay in sequestration
##      onset. This is a discretisation artefact and SHOULD vanish as n_c
##      grows.
##
## They are separated here by window length: (2) acts within a single cycle,
## so it is present even in the first cycle; (1) accumulates over cycles, so
## its contribution grows with the fitted window.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/nc-mechanism.R

suppressPackageStartupMessages({library(plasmofit); library(dplyr); library(tibble)})

CL0   <- 45.012      # truth, the schedule simulation's value
R0    <- 6.5         # within the 5.78-6.78 range of the truth fit
MU    <- 0
T0    <- 1           # log10_total0
DT    <- 0.5
NC_T  <- 96L         # generating n_c

traj <- function(cl, n_c, b_shape, b_off, ts, R = R0, t0 = T0) {
    y0 <- plasmofit:::generate_starts(cl, n_c, b_shape, b_off, t0)
    plasmofit:::mat_exp_series(y0, cl, n_c, R, MU, ts, DT)
}

## ---- channel 2 on its own: the circulating duty cycle ------------------ ##
cat("=== sequestration profile, discretised ===\n")
cat("Cells: mean(y_vals) is the average probability of NOT being sequestered\n",
    "over the cycle, i.e. the circulating duty cycle, at cycle_length =",
    CL0, "h.\n", sep = "")
duty <- bind_rows(lapply(c(96L, 192L, 384L, 1536L), function(nc)
    tibble(n_c = nc, step_h = CL0 / nc,
           duty_cycle = mean(plasmofit:::make_y_vals(CL0, nc)))))
duty$vs_finest <- duty$duty_cycle - duty$duty_cycle[nrow(duty)]
print(as.data.frame(duty), digits = 5)
cat("\nIf the duty cycle falls with n_c, a fit at higher n_c sees a SHORTER\n",
    "circulating window per cycle and must shorten cycle_length to keep the\n",
    "observed fraction of the period that parasites are visible.\n", sep = "")

## ---- the forward map, by window length and synchrony ------------------- ##
best_cl <- function(n_c_fit, b_shape, b_off, tmax, free_nuis) {
    ts <- seq(DT, tmax, by = DT)   # mat_exp_series requires ts > 0
    y_true <- traj(CL0, NC_T, b_shape, b_off, ts)
    lt <- log10(y_true + 1)
    obj <- function(par) {
        cl <- par[1]
        R  <- if (free_nuis) exp(par[2]) else R0
        t0 <- if (free_nuis) par[3] else T0
        yh <- tryCatch(traj(cl, n_c_fit, b_shape, b_off, ts, R = R, t0 = t0),
                       error = function(e) NULL)
        if (is.null(yh) || any(!is.finite(yh))) return(1e12)
        sum((log10(yh + 1) - lt)^2)
    }
    if (free_nuis) {
        o <- optim(c(CL0, log(R0), T0), obj, method = "Nelder-Mead",
                   control = list(reltol = 1e-10, maxit = 2000))
    } else {
        o <- optimize(function(cl) obj(cl), c(35, 50), tol = 1e-6)
        o <- list(par = o$minimum)
    }
    o$par[1]
}

cat("\n=== best-matching cycle_length at a different n_c ===\n")
cat("Cells: trajectory generated at n_c = 96 with cycle_length = ", CL0,
    " h; the\nvalue shown is the cycle_length that best reproduces it at the\n",
    "fitted n_c, in hours, with the shift from ", CL0, " in brackets.\n",
    "`nuisances` says whether R and log10_total0 were re-optimised too.\n\n",
    sep = "")

grid <- expand.grid(b_shape = c(15, 400), tmax = c(48, 96, 144),
                    n_c_fit = c(192L, 384L), free_nuis = c(FALSE, TRUE))
res <- bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    cl <- best_cl(g$n_c_fit, g$b_shape, 0.25, g$tmax, g$free_nuis)
    tibble(b_shape = g$b_shape, window_h = g$tmax, n_c_fit = g$n_c_fit,
           nuisances = if (g$free_nuis) "free" else "fixed",
           best_cl = cl, shift = cl - CL0)
}))
print(as.data.frame(res |> arrange(nuisances, b_shape, n_c_fit, window_h)),
      digits = 4)

cat("\nReading it:\n",
    " - shift NEGATIVE reproduces the real-data direction.\n",
    " - shift roughly CONSTANT across window length -> channel (2), the\n",
    "   within-cycle sequestration discretisation, since that acts in the\n",
    "   first cycle and does not accumulate.\n",
    " - shift GROWING with window length -> channel (1), desynchronisation,\n",
    "   which accumulates over cycles.\n",
    " - larger shift at b_shape 400 than 15 -> the production configuration\n",
    "   is MORE sensitive to n_c, not less.\n", sep = "")
