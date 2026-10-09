# The age-structure decision: chain plus a parameter, or an IPM

Drafted 2026-10-08 while SLURM 30524 was running; all jobs it refers to have
since finished. Thread 15. This is the
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

**That check RAN: SLURM 30527, 2026-10-08, and is reported above.**
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

**This changes the cost calculus** -- or so it seemed on 2026-10-08. The
claim was that the speedup needs no rewrite and no new model, so if cost is
the motive one should port the convolution and leave the IPM to argue only
for decoupling.

> **SUPERSEDED 2026-10-09.** The 38x is against `mat_exp_series`, which the
> fitted model never calls -- it calls `ew_poly_eval`. Measured through
> `grad_log_prob` on the real data the convolution is **5x slower** at
> `n_c` = 96 and 2.8x at 768, and less accurate at high `n_c`. It was added
> to the package and removed again (`bba580f`). The speedup is real only for
> the verification path and the R-side ML screens, which is where it is used.
> The conclusion below is unaffected: the IPM was already rejected on
> identifiability, not on cost.

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

## The profile that decides it is now dense, and running

**SLURM 30576**, submitted 2026-10-08, 14 tasks, ~30–45 min.
`_scripts/nc-profile-fast.R`, read with `_scripts/nc-profile-fast-read.R`.

The convolution makes a whole unit fit cost 12 s instead of 1777 s, so the
profile is no longer limited to three rungs. It now has:

- the **exact chain** at nine integer rungs, 64 to 1024 — the floor question;
- the **gamma-kernel IPM** at a fixed mesh of 192 with sixteen **continuous**
  `n_eff`, 48 to 4096 — the identifiability question;
- a **mesh-convergence check** repeating three rungs at M = 384, so the IPM
  rows can be shown not to be mesh artefacts.

The continuous `n_eff` is the part the discrete ladder structurally could not
provide. In the chain, `n_c` moves the mesh and the rate together, so a profile
over it cannot separate "the dispersion is identified" from "the mesh changed".
Holding the mesh fixed and varying `n_eff` separates them, and that is the only
honest way to ask whether an IPM would return a rate or a boundary value.

Validated before submission: on unit 1 it reproduced 29684's maximised
log-likelihoods at `n_c` 96/192/384 to **1e−6**. The script stops if any chain
rung disagrees with 29684 by more than 0.01 log-likelihood units.

**The reading rule is unchanged** — the one pre-registered in
`_scripts/nc-dispersion-profile-read.R` before any of this existed: turns over
/ saturates / still climbing / flat, judged on whether the 2-log-likelihood
interval is bounded in the fine direction.

## 30576 read, 2026-10-08: no verdict, and what a re-run needs

The job completed cleanly and its forward map checked out — every unit matched
SLURM 29684 to 1e−11 or better — but **both profiles are rougher than the
2-log-likelihood currency the reading rule uses** (2.39 for the chain, 7.17
for the gamma IPM), so the reader returns NO VERDICT and **nothing here
changes the decision above**. Full numbers in `claude/findings.md`.

The thing to carry forward is *why*, because it nearly produced a wrong
answer. The pooled gamma profile dropped 11.10 units at its finest rung, which
reads as a turnover — the branch that would have made dispersion a measured
quantity and been the strongest case for building an IPM. It was **one unit's
optimiser failure**: refitting with eight starts instead of two,
`DSM265|1800` alone gained +11.29.

The mesh-convergence check refit that rung at twice the mesh and moved it by
0.50, which looked like confirmation. It was not: both meshes ran the same
optimiser from the same starts. **A resolution check cannot detect optimiser
error.**

So the identifiability question — does an IPM return a rate or a boundary
value — is still open, and a re-run needs two changes, both affordable now and
neither affordable before the convolution:

1. **Eight starts rather than two.** Measured as necessary, not assumed.
2. **The `b_shape` cap raised or removed.** It binds at 5000 in 13–14 of 14
   units at every rung with `n_eff` ≤ 157, so the whole coarse arm of both
   profiles reports a lower bound rather than a maximum. Note this is a choice
   with content: the production Stan model pins `b_shape` at 400 with
   `max_shape` 1000, so a screen that lets it reach 5000 is already outside the
   production range, and raising it further moves the screen further from the
   model it is meant to inform.

---

## DECISION, 2026-10-09: do not build it

The eight-start profiles settle it. Both **SATURATE**: the data put an upper
bound on dispersion and **no lower bound**, so they do not exclude zero, and
an IPM would return `sigma_d` pressed against zero with an interval touching
it — a boundary estimate rather than a measured rate. That was the branch this
memo identified as arguing *against* building one, and it is the branch that
fired.

| profile | best | 2-unit interval | top-rung gain | unimodality violation |
|---|---|---|---|---|
| exact chain, `n_c` 64–1024 | 1024 | {512, 1024} | +1.30 | 0.00 |
| gamma IPM, mesh fixed at 192 | `n_eff` 2264 | {1683, 4096} | +0.69 | 0.69 |

**Every argument in this memo has now resolved against the rewrite.**

- **Cost** went first: the exact chain is a convolution, 38× faster per
  trajectory and 148× on a unit fit, with no rewrite and no change of model.
- **The decay law** could not distinguish √ from linear, so it never supported
  a kernel choice.
- **Identifiability**, the last prop, has come back negative.
- **Option A (non-uniform stage rates) is also dead**, and for the same
  reason: it can only *add* dispersion above the Erlang floor, and the data
  want less, not more.

**Option C is the answer.** Pick `n_c` by elpd, report `cycle_length` with the
ladder, and state the sensitivity — which is now small: 40.93 to 41.15 h
across the 2-log-likelihood interval once `n_c` ≥ 384, against 3.20 h over the
original ladder.

**And the underlying reason is worth more than the decision.** The data
constrain the **total** stage spread at the end of the window, not how it
splits between initial synchrony and accumulated desynchronisation: the
accumulated term varies 2.0–4.4× across rungs while the total varies 12.7–23.8%,
and `b_shape` falls monotonically with the rung in 14 of 14 units. No
reparameterisation can recover a quantity the data do not separately identify.
A rewrite would have moved the problem, not solved it.
