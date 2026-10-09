// NOTE, 2026-10-09: needs a plasmofit built from commit 05c9c1f;
// conv_series() was removed from the package in bba580f.
// Gradient and log-density check for conv_series(), against mat_exp_series().
//
// Parsing proves nothing about reverse-mode autodiff through a complex FFT, so
// this compiles both paths behind a data switch and lets the driver compare
// log_prob and grad_log_prob at identical parameter values. mat_exp_series()
// is the reference: its gradients come from Stan's own matrix_exp, which is
// well tested, so agreement between the two is a real check on the FFT adjoint
// rather than a self-consistency check.
//
// Driven by _scripts/conv-series-grad.R.
functions {
#include functions/plasmofit.stan
}
data {
    int<lower=1> n_ts;
    array[n_ts] real ts;
    vector<lower=0>[n_ts] y_obs;
    int<lower=1> n_c;
    int<lower=1> M;                  // FFT length, from conv_min_M()
    real<lower=0> dt;                // step for mat_exp_series
    int<lower=0, upper=1> use_conv;  // 1 = conv_series, 0 = mat_exp_series
}
parameters {
    real<lower=35, upper=50> cycle_length;
    real<lower=2.5, upper=4000> b_shape;
    real<lower=0, upper=1> b_offset;
    real<lower=2, upper=9> log10_total0;
    real<lower=1.5, upper=100> R;
    real<lower=0> sigma;
}
model {
    vector[2 * n_c] y0 = generate_starts(cycle_length, n_c, b_shape,
                                         b_offset, log10_total0);
    vector[n_ts] y_hat;
    if (use_conv == 1) {
        y_hat = conv_series(y0, cycle_length, n_c, R, 0.0, ts, M);
    } else {
        y_hat = mat_exp_series(y0, cycle_length, n_c, R, 0.0, ts, dt);
    }
    target += normal_lpdf(log10(y_obs + 1) | log10(y_hat + 1), sigma);
    target += lognormal_lpdf(sigma | log(0.3), 0.5);
    target += normal_lpdf(cycle_length | 45, 5);
}
