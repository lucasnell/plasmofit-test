# The age-structure decision: chain plus a parameter, or an IPM

Drafted 2026-10-08, while SLURM 30524 is running. Thread 15. This is the
design document behind the `TODO.md` item "Give the desynchronisation rate its
own parameter"; `claude/threads.md` thread 15 holds the running record and
`claude/findings.md` holds every number quoted here.

The case is written before 30524 reads out, deliberately, so that the reading
rule is not reverse-engineered from the answer.

## The defect, stated precisely

The Erlang chain has `n_c` sequential compartments with common rate
`lambda = n_c / cycle_length`. Transit through the chain is therefore
Erlang(`n_c`, `lambda`):

- mean = `cycle_length`
- variance = `cycle_length^2 / n_c`
- **coefficient of variation = `1 / sqrt(n_c)`**

and the stage distribution's sd after `k` cycles is
`cycle_length * sqrt(k / n_c)`.

So one integer sets two unrelated things:

1. **The numerical mesh.** The resolution at which developmental age is
   represented, and in particular the resolution at which the sequestration
   logistic is evaluated — `make_log_y_vals` puts the circulating fraction on
   a logistic in absolute age with `p3 = 18.5802` h, sampled at
   `k * cycle_length / n_c`.
2. **The biological desynchronisation rate.** How fast an initially
   synchronous cohort spreads out.

There is no reason these should be the same number, and raising `n_c` to
improve (1) suppresses (2).

## What is established

| | evidence |
|---|---|
| `n_c` = 96 is wrong | 71.8 elpd (se 13.2) behind 192 under Bayes; **+72.4 log-likelihood units in 14 of 14 units** under a priors-free ML profile |
| Not a numerical artefact | `max_rel_diff` is 8e−13 at `n_c` = 96 and **grows** with `n_c`, so better arithmetic is not what the higher rungs are buying |
| It costs real hours | Misspecifying `n_c` moves `cycle_length` by 1.0–1.45 h in either direction (2×2 simulation, 29637); the real ladder spans 3.20 h |
| The observable period is not the parameter | A one-sided sequestration window (~40% duty cycle) convolved with a broadening right-skewed age distribution shifts the apparent period; ~0.75 h short at `n_c` = 96 |
| The decay **law** is not identified | At the predicted sigma = 0.033 a linear-in-cycles component buys at most +0.396 and typically +0.0000 log-likelihood units over 1130 observations (29684) |

**The pattern across the last two results is the key to this decision.** The
data pin the *amount* of dispersion sharply — that is what a 71.8 elpd gap
between `n_c` 96 and 192 is, since both rungs obey the same √ law and differ
only in magnitude — and they cannot resolve its *functional form* at all.

So the deliverable of any fix is a **free magnitude**. Expecting the data to
teach us the form is not supported, and designing for it is designing for
something 29684 showed is invisible here.

## The floor, which is what the whole decision turns on

For a continuous-time Markov chain on `n` sequential states with rates
`lambda_1 .. lambda_n`, transit time is the sum of independent exponentials:

- mean = `sum(1 / lambda_i)`
- variance = `sum(1 / lambda_i^2)`

Minimising `sum(x_i^2)` subject to `sum(x_i) = mean` puts every `x_i` equal, so
**the Erlang is the minimum-variance member of the family** and
`CV >= 1 / sqrt(n)`. Making the rates unequal only ever *adds* variance, up to
`CV -> 1` when one slow stage dominates.

This is a property of the whole chain family, not of the equal-rate special
case. A deterministic delay — the obvious way to shrink variance — is Erlang
with infinitely many stages, so it is not available at finite `n`.

**Therefore the one thing no chain can do, at any affordable `n_c`, is produce
dispersion below `1 / sqrt(n_c)`.** Everything above the floor is reachable
with one extra parameter. That is exactly the question SLURM 30524 asks.

## Option A — keep the chain, free the stage rates

Replace the common `lambda` with a rate profile carrying one shape parameter,
holding `sum(1 / lambda_i) = cycle_length`. A two-level profile (a fraction of
stages fast, the rest slow) or a monotone profile both work; the transit
distribution becomes hypoexponential rather than Erlang.

**What it changes:** the generator's diagonal. Nothing else. The generator
stays bidiagonal, `mat_exp_series` is unchanged, and the `matrix_exp`
cross-check still applies — it needs revalidating at non-uniform rates, but
the instrument survives.

**What it buys:** `n_c` can then be chosen for the mesh alone, and dispersion
dialled anywhere from the floor upward. Since mesh accuracy is already
excellent at `n_c` = 96 (`max_rel_diff` 8e−13), the binding mesh requirement
is the sequestration-window resolution, which is an O(1/`n_c`) effect worth
0.284 h between 96 and infinity. Pick `n_c` = 384 for that, giving a floor of
CV = 0.051, and the whole range above it is free.

**What it cannot do:** go below the floor.

**Compatibility:** Erlang is the boundary case, so **every existing fit remains
nested and comparable**. This matters more than it sounds — the project's
evidence is a web of paired comparisons against earlier fits.

**Cost:** one parameter, a change to how the generator is built in Stan, a
revalidation of the cross-check, and a refit of the ladder.

**Risk:** if the data want the floor, the new parameter sits at its boundary.
That is a **result, not a failure** — it would say the data want the least
dispersion the chain can give at that mesh — but it would also mean the
parameter is inert and the reported cycle length still depends on `n_c`
through the floor.

## Option B — replace the chain with transport plus dispersion (IPM)

Carry a density over continuous developmental age on `[0, cycle_length)`. Each
step advects by `dt`, convolves with a dispersion kernel of width `sigma_d`,
and multiplies by `R` at the cycle boundary. The mesh `M` is chosen purely for
quadrature accuracy and `sigma_d` is a free parameter independent of it.

**What it buys:** decoupling in **both** directions. `sigma_d -> 0` is
deterministic transit at any mesh, which no chain can reach.

**What it does not buy:** a different decay law. With a Gaussian kernel the sd
after `k` cycles is `sigma_d * sqrt(k)` — the same √ law. This is a change of
parameterisation, not of mechanism. 29684 says that is acceptable, because the
data cannot distinguish the laws anyway, but it must not be sold as testing
the form.

**The implementation hazard, which is the serious one.** An advection scheme
on a mesh introduces **numerical diffusion** — artificial spreading that is a
function of the mesh, which is precisely the defect being removed, re-entering
through the back door. Any IPM here has to use a scheme where that is
controlled and measured, not assumed small. The natural choice on a periodic
age domain is spectral: advection is a phase rotation and Gaussian convolution
is a multiplication, both exactly diagonal in Fourier space, so neither step
adds mesh-dependent spread. That is a clean design, and it should be the
default proposal rather than a finite-volume scheme.

**Cost, and it is the real objection.** Discards the validated Erlang-window
series and its `matrix_exp` cross-check (good to 1e−12), requires revalidating
the entire numerical path, rewrites `generate_starts` and `mat_exp_series` and
therefore the schedule-simulation machinery, needs new priors, and **makes
every existing fit incomparable**. This is a project-scale change, not a
parameter change.

**What it buys scientifically, which is the real argument for it.** The
desynchronisation rate becomes a reported estimate with an interval instead of
a fixed assumption that moves `cycle_length` by 1.0–1.45 h. "Parasites
desynchronise at rate X per cycle, 95% CI [a, b]" is a result; "we set
`n_c` = 192" is an assumption the reader has to take on trust. That is a
manuscript-level argument, not a likelihood argument, and it should be made as
such rather than dressed up as better fit.

## Option C — change nothing structural

Pick `n_c` by elpd, report `cycle_length` conditional on it, and publish the
ladder as a sensitivity analysis.

**Cost:** zero. **What it concedes:** the headline number carries an
assumption-driven uncertainty of 1.0–1.45 h, comparable to the bias the
project has spent months attributing. The honest form of C is to report the
ladder in the main text rather than a supplement.

## How 30524 decides between them

Pre-registered in `_scripts/decay-law-768.sh` and enforced by
`_scripts/nc-768-read.R`, which prints the verdict itself.

| summed gain 384→768 | units positive | verdict | implication |
|---|---|---|---|
| < +3 | < 9/14 | plateau | the preferred dispersion is **at or above** the floor → **A** suffices, and is far cheaper |
| ≥ +8 | ≥ 10/14 | still climbing | the data want dispersion **below** the floor → only **B** can deliver it |
| anything else | | ambiguous | decide on biology; **C** is the honest default |

Geometric decay of the earlier gains (+72.4 then +10.4, ratio 0.14) predicts
+1.5, which falls in the plateau branch. That is a prediction, not a
preference.

## Two things that hold whichever way it reads

**The form will not be learned from these data.** 29684 settled that. Whatever
is built, the kernel's shape is an assumption, and the manuscript has to say
so. Only the magnitude is estimable.

**Identifiability does not vanish under B.** `b_shape` sets the spread at
t = 0 and `sigma_d` sets its growth; with observations from 1.6 to 4.8 cycles
these are an intercept and a slope on the same curve, so they separate in
principle. But 29684's lack of power is a direct warning that the slope is
weakly determined over this window, and `sigma_d` would likely need a prior
doing real work. **Check that before building**, with the same cheap ML
machinery: profile `sigma_d` with `b_shape` free and see whether the profile
has curvature. If it does not, B delivers a parameter the data cannot
estimate, and the honest outcome is C with a clear caveat.
