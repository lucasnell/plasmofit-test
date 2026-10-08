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

**That check is running: SLURM 30527**, submitted 2026-10-08.
`_scripts/decay-law-profile.sh` adds `n_c` = 128, 256, 512 so the profile has
seven rungs spanning a transit CV of 0.102 down to 0.036, and
`_scripts/nc-dispersion-profile-read.R` reads it. It also records the fitted
`b_shape` at every rung, which the earlier output dropped, so the
`b_shape`/dispersion trade-off — the actual confound — can be read directly.

The rule was refined before any new rung existed, because a dry run on the
three rungs already on disk showed the first draft asked the wrong question.
It tested for an **interior maximum**, and a profile that saturates toward an
asymptote never has one. What matters is whether the 2-log-likelihood interval
is **bounded in the fine direction** — whether the data exclude zero
dispersion.

**Note which way saturation cuts.** If the profile saturates, the data put an
upper bound on dispersion and **no lower bound**, so an IPM would return
`sigma_d` pressed against zero with an interval touching it — a boundary
estimate, not a rate. That argues **against** building one, the opposite of
what the structural argument alone suggests. On the three rungs already on
disk the reader says *still climbing*: the top rung still buys +10.36
log-likelihood units, which is why both 30524 and 30527 exist.

## What the prototype found, 2026-10-08

`_scripts/ipm-prototype.R`, output `_data/ipm-prototype-2026-10-08.txt`.
Three results, two of which correct this memo above.

### The exact chain is a convolution, and it is 38x faster

Work in **absolute** developmental age rather than age modulo the cycle. Then
transport is a pure-birth process — a parasite's absolute stage after time `t`
is its starting stage plus Poisson(`lambda * t`) — with no wrap and no
boundary condition. Growth becomes a weight, `R^d` for `d = floor((j-1)/n_c)`
divisions, and sequestration becomes a weight too, because circulating status
is a function of position within the cycle and **resets at division** (the ODE
sends both `N[n_c]` and `S[n_c]` into `N[1]`).

So the whole trajectory is a convolution evaluated once per observation time.
It reproduces `mat_exp_series` to **1e-12**, the same order as the package's
own `matrix_exp` cross-check, and runs **38x faster at `n_c` = 384** (1.179 s
against 0.031 s). It is a reimplementation from the ODE's structure, not a
refactor, so the agreement is an independent check of both.

**This changes the cost calculus.** The speedup needs no rewrite and no new
model: it is the same model, agreeing to 1e-12. If cost is the motive, port
the convolution. The IPM is then only for decoupling, which is the honest way
to argue it.

Two subtleties the validation caught, both of which a fresh implementation
would plausibly get wrong in silence:

1. The chain applies the sequestration hazard with a **one-stage lag** — the
   transition from stage `k-1` to `k` uses `q[k-1]` — so the circulating
   fraction is `G[k-1]`, not `y[k]`. This is the O(1/`n_c`) delay the Stan
   source documents. Using `y[k]` costs 12% at the first observation.
2. Parasites in their **first cycle** have not been reset, so their weight
   depends on where they started. It is separable, so it costs one extra
   convolution rather than breaking the method.

### Correction: an IPM is not a reparameterisation of the chain

The memo says above that a Gaussian kernel gives the same sqrt law, so an IPM
"is a change of parameterisation, not of mechanism". **That is wrong.** The
chain's stage at time `t` is exactly **Poisson — a lattice distribution** —
and any continuous kernel differs from it in the tails.

Measured at matched mesh and matched variance, in log10 units, which is the
scale the likelihood works on (the residual sd of these fits is about 0.48):

| `n_c` | gaussian | gamma |
|---|---|---|
| 96 | 0.0135 | 0.0120 |
| 192 | 0.0331 | 0.0310 |
| 384 | 0.0623 | 0.0595 |

**The gap grows with `n_c`, and skew is not the cause.** A gamma kernel
preserves the chain's right skew exactly — `Gamma(shape = n_eff*t/cl, scale =
cl/n_eff)` is the chain with *continuous* `n_c` — and it lands within 0.003 of
the Gaussian. The cause is the troughs: the observable spans four orders of
magnitude across the cycle because sequestration hides ~60% of it, and at a
trough the value is set by the **tail** of the age distribution, which is
exactly where a lattice Poisson and a continuous kernel part company. At
`n_c` = 384 the entire 0.0595 sits at **one** time, `t` = 72 h, the deepest
trough; the other six agree to 0.006 or better.

Troughs are also where **low-end censoring** is an open question in `TODO.md`,
so this is the one place where the kernel choice and a known data issue
coincide. Any IPM fit has to be compared with the chain at the troughs
specifically, not on an average.

### Correction: the numerical-diffusion hazard does not arise

The memo recommends a spectral scheme because advection on a mesh adds
artificial spreading. **There is no time stepping at all**, so there is
nothing to accumulate: the kernel is applied analytically, once per
observation time. The spectral recommendation is superseded by something
simpler and strictly better.

### The decoupling works, and here is the number

Max |log10 difference| from the chain at `n_c` = 384, which is the dispersion
the data prefer:

| forward map | difference |
|---|---|
| chain at `n_c` = 96 — mesh and rate locked together | **1.589** |
| IPM gamma, mesh M = 96, `n_eff` = 384 | 0.064 |
| IPM gamma, mesh M = 192, `n_eff` = 384 | 0.061 |

A factor of 25. A coarse mesh with the fine dispersion dialled in is a close
match to the fine chain, while the coarse chain is not — mesh and rate have
been separated, which is the thing the whole design argument rests on.

### What this does to the three options

**Option A** (non-uniform stage rates) is unchanged and still cheapest.

**Option B** (IPM) is cheaper to build than this memo assumed, because the
convolution machinery is already written and validated, and `n_eff` in the
gamma kernel is a one-line continuous parameter. But it is a **different
model**, not a reparameterisation, so every fit would have to be redone and
compared at the troughs.

**Option C** (change nothing) gains a new argument: port the Poisson
convolution, keep the model exactly as it is, and spend the 38x on more
replicates and a finer `n_c` ladder rather than on a rewrite.
