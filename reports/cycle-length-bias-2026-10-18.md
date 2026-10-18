# Why the model's cycle-length estimates are biased: current understanding

*plasmofit, 18 October 2026. Draft for sharing. Every number below names the
script and saved output it comes from (repository
`github.com/lucasnell/plasmofit-test`; scripts in `_scripts/`, outputs in
`_data/`, which is not tracked). Each claim carries a strength label:
**Established** (direct test, reproduced), **Supported** (direct test, thin
replication or a known caveat), or **Open** (not yet tested at the settings
that matter).*

## Summary

1. **The bias exists.** On data simulated from known parameters at the real
   observation schedules, the model returns a cycle length that is too long:
   +1.97 h on average over 3 simulated datasets (range +1.14 to +2.58 h,
   truth 45.01 h), with 96 developmental compartments used both to simulate
   and to fit. It survives corrected priors (+1.63 to +1.83 h) and is still
   +1.55 h at 192 compartments. **Established.**
2. **Most of it is the cost of estimating the nuisance parameters, and among
   them only the initial-synchrony parameter matters.** With every nuisance
   parameter held at its true value, the likelihood alone peaks +0.41 h high.
   Fixing initial synchrony at the truth removes about 0.9–1.0 h; the
   location of the cycle-length prior (centred at 48 h) accounts for about
   0.44 h. The other nuisance parameters contribute nothing. **Supported**
   (the decisive fixed-synchrony arm has one replicate that passes the
   convergence gate; a weaker version of the same intervention, on six
   replicates, agrees in direction).
3. **The underlying reason is an identification gap.** Parasites are first
   observed 1.6 cycles after inoculation, so the data constrain only the
   total spread of developmental ages during the observation window, not how
   that spread divides between synchrony at inoculation and
   desynchronisation accumulated since. Initial synchrony is therefore not
   identified, and the fit trades oscillation amplitude against period: an
   error in synchrony is absorbed partly by the cycle length. **Established**
   for the identification gap; **Supported** for the trade-off as the route
   to cycle length.
4. **A second, separate source: the number of compartments is a biological
   assumption.** It fixes the desynchronisation rate. Because only about the
   first 40% of each cycle is visible before parasites sequester, the
   observable period of the oscillation differs from the cycle-length
   parameter by an amount that depends on the compartment count (0.75 h
   short at 96, within 0.25 h at 384). Fitting with a compartment count wrong by
   a factor of two moves the estimate by 1.0–1.5 h. On the real data the
   estimate fell from 44.24 h to 41.04 h as the count rose from 96 to 384.
   **Established** for the mechanism; **Supported** for its size on real
   data.
5. **Ruled out:** solver mismatch, the posterior summary statistic, the
   prior's width, the parameter bounds, the initial-density prior, the
   observation-error scale, and two phase-based explanations (period–phase
   compensation and excess phase-space volume at periods commensurate with
   the 12 h sampling grid).
6. **Not yet known:** the size of the bias at the settings now used for the
   real data (384 compartments, synchrony estimated). A simulation study to
   measure it is running; its partial read has no readable cell. **The
   current real-data estimates, 41.3 h and 41.8 h (mean over 13 trials), are
   uncorrected and should not be reported as final.**

## The model, and the terms used below

A hierarchical Bayesian ODE model, fitted in Stan through the `plasmofit` R
package, of parasite density (infected red blood cells per mL) in 13
volunteer infection trials: 1130 observations on 177 series. Observations run
from 72 h to 216 h after inoculation (1.6 to 4.8 cycles), 4–8 per series,
mostly 12 h apart.

| term | meaning |
|---|---|
| cycle length (`cycle_length`) | duration of the asexual blood-stage cycle, in hours; one per trial |
| compartment count (`n_c`) | parasites pass through `n_c` sequential exponential compartments per cycle, so the time to complete a cycle is Erlang-distributed with coefficient of variation 1/√`n_c`. After *k* cycles the spread of developmental age is √(*k*/`n_c`) cycles. `n_c` therefore sets the **desynchronisation rate**, and the model has no separate parameter for it |
| initial synchrony (`b_shape`) | concentration of the developmental-age distribution at inoculation; larger means more tightly synchronised |
| initial phase (`b_offset`) | median developmental age at inoculation |
| sequestration | parasites leave circulation from 18.58 h of development, so only the first ~40% of a ~45 h cycle is observable |
| `R`, `log10_total0`, `sd_iRBC` | multiplication per cycle; log10 initial density; observation-error SD on the log10 scale |
| nuisance parameters | the five above (`b_shape`, `b_offset`, `R`, `log10_total0`, `sd_iRBC`): everything except cycle length |
| convergence gate | a fit counts only if the maximum split R-hat over all parameters is below 1.05 |
| elpd | expected log pointwise predictive density, by PSIS leave-one-out over the 1130 observations; higher is better |

## How bias is measured

Three instruments, which separate different causes:

1. **Simulation and refit.** Simulate data from known parameters at the real
   observation times and group structure, add observation noise, fit the
   model, and compare. **Bias** throughout is: for each trial, posterior mean
   cycle length minus the true cycle length; then the mean over the 13
   trials; in hours; positive means too long. Interventions (fix a parameter,
   move a prior) are refitted to the **same** simulated datasets, so their
   effects are paired differences. Pairing matters: different noise draws of
   the same design move the bias by about 1 h, more than most of the effects
   measured.
2. **Prior-free profile likelihood.** Maximum likelihood, with chosen
   parameters held at the truth, to see where the likelihood itself peaks.
3. **Deterministic forward-map checks.** Noiseless trajectories, no
   estimation, to isolate what is built into the model's arithmetic.

## 1. The model overestimates cycle length in simulation

Simulated from the posterior means of a fit that shares one cycle length
across all trials (45.01 h), at `n_c` = 96, under the priors then in use;
fitted with the same model.

*Cells: bias in hours (mean over 13 trials of posterior mean minus truth);
one simulated dataset per row. Scripts `wockner-schedule-sim.R` and
`wockner-schedule-sim-analyze.R`; output `schedsim-analyze-2026-10-06.txt`.*

| simulated dataset | bias |
|---|---|
| replicate 1 | +2.58 |
| replicate 2 | +1.14 |
| replicate 3 | +2.19 |
| **mean** | **+1.97** |

It is robust to the choices that might have manufactured it:

- **Not a solver mismatch.** The simulator and the likelihood use different
  numerical routines; they agree to a maximum relative difference of 5e−13.
- **Not the truth construction.** Simulating from single posterior draws
  instead of the posterior mean vector gives +1.68 and +0.86 h, inside the
  spread above (`schedsim-analyze.log`).
- **Not the posterior summary.** Mean, median, and mode differ by at most
  0.09 h on a bias of ~2 h (`summary-stat-check.R`,
  `summary-stat-check.log`).
- **Not the priors since corrected.** With the initial-density and
  synchrony priors both widened, and data simulated from a fit under those
  priors, the bias is +1.69 to +2.05 h (mean +1.83, 3 datasets); from the
  original truth under the same priors, +1.63 h (5 datasets)
  (`schedsim-truthfit-check.R`, `schedsim-truthfit-check.log`).
- **Not specific to `n_c` = 96.** Simulated and fitted at 192, the bias is
  +1.55 h (3 datasets; `nc-2x2-read.R`, `nc-2x2-2026-10-07.txt`).

## 2. Where it comes from: estimating initial synchrony

### The likelihood carries about a fifth of it

With every nuisance parameter at its true value and only a common phase shift
free, the profile likelihood over cycle length peaks close to the truth.

*Cells: hours. "Likelihood peak" is where the prior-free profile likelihood
peaks, minus the 45.01 h truth (grid 0.25 h). "Full fit" is the bias of the
full Bayesian fit to the same simulated dataset. Paired within a row.
`phase-volume-test.R`, `phase-volume-test.log`.*

| simulated dataset | likelihood peak | full fit |
|---|---|---|
| replicate 1 | +0.74 | +2.58 |
| replicate 2 | −0.01 | +1.14 |
| replicate 3 | +0.49 | +2.19 |
| **mean** | **+0.41** | **+1.97** |

The likelihood's share varies with the noise draw and tracks the full bias
(largest in replicate 1, smallest in replicate 2). The remaining ~1.56 h
appears only when the nuisance parameters are estimated.

### Which nuisance: only initial synchrony

*Cells: paired change in bias, in hours: the intervention's bias minus the
baseline fit's bias on the same simulated dataset, averaged over n pairs in
which both fits pass the convergence gate. Negative means the intervention
reduces the bias. Baseline bias +1.97 h. The rows are single interventions
against a common baseline and **do not add**. `wockner-schedule-sim-analyze.R`,
outputs `schedsim-analyze-2026-10-05.txt` and `-2026-10-06.txt`.*

| intervention | n | change | range |
|---|---|---|---|
| synchrony (`b_shape`) fixed at its true value | 1 | **−0.92** | — |
| synchrony prior widened (log-scale SD 0.5 → 1.5) | 6 | **−0.51** | −0.87 to −0.23 |
| all nuisance priors widened | 6 | −0.46 | −1.08 to −0.07 |
| cycle-length prior centred on the truth | 3 | −0.44 | −0.51 to −0.29 |
| cycle-length prior centred 2.9 h below the truth | 3 | −0.64 | −0.76 to −0.51 |
| hierarchy on cycle length removed | 2 | −0.26 | −0.29 to −0.23 |
| initial-density prior widened | 6 | +0.15 | −0.16 to +0.43 |
| observation-error SD fixed at its true value | 2 | +0.10 | +0.03 to +0.18 |
| cycle-length bounds moved, width kept | 4 | +0.03 | — |

- **Fixing synchrony is the largest lever.** Only one of three fixed-synchrony
  fits passes the gate, which is why n = 1. The other two fail narrowly
  (R-hat 1.055 and 1.052 after reseeding) and give −1.26 and −0.96 h; all
  three together average −1.04 h. Fixing a parameter made sampling harder,
  not easier, which is itself a sign that synchrony was absorbing something.
- **The well-converged corroboration is the widened synchrony prior**: six
  pairs, all negative. Loosening the prior buys about half of what fixing the
  value buys, a dose–response in the expected direction.
- **No other nuisance moves cycle length.** Widening all nuisance priors does
  no better than widening the synchrony prior alone; the initial-density
  prior and the error SD have no effect.
- **Prior location matters more than prior width.** The cycle-length prior
  is centred at 48 h (motivated by reports of 24 h-multiple cycles under
  circadian coordination). Widening it barely moves the estimate because its
  median stays at 48 h; moving its centre to the truth removes 0.44 h, and
  the effect saturates. A correctly placed prior still leaves +1.53 h.

Taken together: likelihood +0.41 h, synchrony estimation ~0.92 h, and prior
location ~0.44 h sum to 1.77 of 1.97 h. Because the last two overlap (both
remove uncertainty), this says the pieces are the right size to close the
budget, not that it decomposes exactly.

### Why synchrony is not identified

**Established.** The first observation is 1.6 cycles after inoculation. By
then the developmental-age spread is the sum (in quadrature) of the spread at
inoculation, set by `b_shape`, and the spread accumulated since, set by
`n_c`. The data see the total and not the split.

The evidence, from a prior-free likelihood profile over the desynchronisation
rate (`nc-profile-fast.R`, `nc-profile-fast-2026-10-09.txt`): across settings
where the accumulated spread varies 2.0-fold (compartment model) or 4.4-fold
(a continuous-age alternative), the fitted **total** spread at the last
observation varies by only 12.7% and 23.8% of its mean. The fitted `b_shape`
falls as `n_c` rises in all 14 trial-by-inoculum units (Spearman
correlation −1.00, median over units): the fit holds the total fixed and
trades one component against the other.

The same gap appears in simulation as attenuation. *Cells: `b_shape` mean
over the 14 trial-by-inoculum units, true value and posterior mean; one
simulated dataset per value in the range. `schedsim-truthfit-check.log`.*

| true `b_shape` | recovered, range over datasets |
|---|---|
| 14.9 (5 datasets) | 19.5 to 27.9 |
| 65.0 (3 datasets) | 32.2 to 37.4 |

A 4.4-fold change in the truth produces a 1.5-fold change in the estimate:
both are pulled toward ~25. On the real data the posterior for `b_shape` is
barely narrower than its prior (coefficient of variation 0.43 against 0.53),
and widening the prior moves the posterior mean from 14.9 to 65.4
(`panel-bshape-spread.log`).

### How synchrony error reaches cycle length

**Supported.** The model can reproduce an observed oscillation with high
amplitude and a short period, or with lower amplitude and a longer one.
Synchrony sets amplitude; if it is misjudged, the period adjusts to
compensate. On the real data, with `n_c` fixed, widening the synchrony prior
moved the mean cycle length over 13 trials by −0.47 h, and by −1.29 h once
the initial-density prior was also corrected (`wockner-prior-panel.R`,
`prior-panel.log`). These are shifts on real data, where there is no truth,
so they show sensitivity, not bias.

## 3. A second source: the compartment count sets the desynchronisation rate

### The observable period is not the cycle-length parameter

**Established.** Deterministic, noiseless trajectories with no estimation.

*Cells: the third peak-to-peak interval of the detrended log10 trajectory
minus the cycle length that generated it, in hours; negative means the
visible oscillation runs short of the parameter. Peaks refined below the
0.25 h grid. `nc-period-check.R`, `nc-period-check-2026-10-07.txt`.*

| `b_shape` | `n_c` = 96 | 192 | 384 |
|---|---|---|---|
| 15 | −0.77 | −0.28 | −0.05 |
| 400 | −0.74 | −0.14 | +0.24 |

Mechanism: only parasites younger than 18.58 h circulate, so the visibility
window is one-sided. As the age distribution broadens and becomes
right-skewed at a rate set by `n_c`, the centre of the visible subpopulation
drifts relative to the true cycle boundary, and the drift reads as a change
in period. The intervals shorten from cycle to cycle (44.62, 44.35, 44.24 h at
`n_c` 96, `b_shape` 15), as this account requires. A fit with too few
compartments must therefore lengthen the cycle-length parameter to match an
observed period. A competing channel (the coarseness of the sequestration
grid at low `n_c`) was checked and has the wrong sign and an eighth of the
size (`nc-mechanism.R`, `nc-mechanism-2026-10-07.txt`).

Matching the whole trajectory by least squares rather than peak intervals,
the shift needed over a 96 h window is −1.0 to −2.4 h depending on `b_shape`
and the target `n_c` (same output).

### Misspecifying the compartment count moves the estimate both ways

*Cells: bias in hours (mean over 13 trials of posterior mean minus the 45.01 h
truth), averaged over simulated datasets that pass the gate (n in brackets).
The diagonal is correctly specified. `nc-2x2-read.R`,
`nc-2x2-2026-10-07.txt`.*

| | fitted at `n_c` 96 | fitted at `n_c` 192 |
|---|---|---|
| **simulated at 96** | +1.97 (3) | +0.93 (3) |
| **simulated at 192** | +3.00 (2) | +1.55 (3) |

Fitting with too few compartments adds +1.45 h over the correctly specified
fit; too many removes 1.04 h. So an assumption about the desynchronisation
rate that is wrong by a factor of two moves the estimate by about an hour, in
the direction the mechanism predicts.

### What this did to the real-data estimates

*Cells: posterior means, averaged over the 13 trials (cycle length) or the 14
trial-by-inoculum units (`b_shape`); elpd summed over 1130 observations.
`b_shape` estimated under the widened prior in all three fits.
`nc-ladder-read.R`, `nc-ladder-2026-10-06.txt`.*

| `n_c` | `b_shape` | cycle length (h) | elpd | max R-hat |
|---|---|---|---|---|
| 96 | 64.45 | 44.24 | −926.2 | 1.019 |
| 192 | 61.57 | 42.07 | −857.6 | 1.032 |
| 384 | 25.94 | 41.04 | −854.4 | 1.055 |

*Paired elpd differences, higher count minus lower, summed over the 1130
observations, with the standard error of the paired difference
(`nc-ladder-loo.R`, `nc-ladder-loo-2026-10-18.txt`):* 192 over 96, +68.6
(9.7); 384 over 96, +71.8 (13.2); 384 over 192, +3.2 (5.0). So 96 is clearly
worse, and 192 and 384 are indistinguishable. This is not numerical error:
the model's own cross-check puts relative error at 8e−13 at 96, and it grows
slightly with `n_c` (`nc-numerical-check-2026-10-06.txt`). It is the biology
the model assumes, and it holds without priors: a prior-free profile
likelihood favours 192 over 96 by 72.4 log-likelihood units summed over the
14 trial-by-inoculum units, with all 14 agreeing (`decay-law-read.R`,
`decay-law-read-2026-10-08.txt`). The 384 fit sits marginally above the
gate; the comparison that carries the conclusion is 96 against 192.

Above 384 the estimate stabilises: across the range of compartment counts
within 2 log-likelihood units of the best in the prior-free profile, cycle
length spans 40.93 to 41.15 h (`nc-profile-fast-2026-10-09.txt`).

Roughly half to two thirds of the real-data fall from 96 to 192 or 384 compartments
is reproduced by the noiseless forward-map effect alone. The remainder is not
attributed.

## 4. One phenomenon, seen twice

Sections 2 and 3 are the same trade-off approached from two sides. The data
constrain the total developmental-age spread over the window. How that spread
is attributed, to synchrony at inoculation or to desynchronisation since,
changes the amplitude the model can produce at each time and how far the
visible oscillation's period departs from the cycle-length parameter. The
cycle length absorbs the mismatch. With `n_c` correctly specified, estimating
the unidentified synchrony costs about an hour; with `n_c` misspecified, the
fixed desynchronisation rate forces a further error of similar size. Both
push in the same direction at the settings this project began with, which is
why the early estimates (~45 h) were high.

This is the same interaction Greischar and Childs (2023, *Trends in
Parasitology* 39(8), doi 10.1016/j.pt.2023.05.006) describe for the parasite
multiplication rate: some developmental ages are easier to sample than
others, and the age distribution changes over the infection. They report the
largest errors in synchronised infections and when the median parasite age is
offset from the sampling times by about 12 h. Here 87.8% of within-series
sampling intervals are 12 h. Their result concerns multiplication rate and a
different estimator, so it motivates a prediction rather than transferring
directly: the cycle-length bias should depend on the true initial phase
relative to the sampling grid and should be larger at high synchrony.
**Open**: every simulation so far shares one true phase.

## 5. Explanations tested and ruled out

*Cells: paired change in bias in hours against the baseline fit, as in the
table in section 2, or the stated test result.*

| explanation | test | result |
|---|---|---|
| simulator and likelihood use different solvers | maximum relative difference over 150 draws | 5e−13 |
| skewed posterior, wrong summary | mean vs median vs mode | differ by ≤ 0.09 h |
| cycle-length prior width | prior variance quadrupled | −0.15 h |
| bounds on cycle length | bounds moved, width kept (n = 4) | +0.03 h |
| initial-density prior, or its trade-off with `R` | prior widened (n = 6) | +0.15 h, sign inconsistent |
| observation-error scale | fixed at truth (n = 2) | +0.10 h, wrong sign |
| period and initial phase compensate along a ridge | posterior correlation; phase decomposition | \|r\| ≤ 0.16; no ridge |
| more phase-space volume at 48 h (commensurate with the 12 h grid) | profile vs phase-integrated likelihood | both peak at the same place; volume at 48 h is slightly lower |
| sampling schedules alias the likelihood | per-trial maximum-likelihood estimates | mean bias −0.20 h, indistinguishable from zero |

Sources: `wockner-schedule-sim-analyze.R`, `summary-stat-check.R`,
`schedsim-phase-period.R`, `phase-volume-test.R`,
`wockner-schedule-bias-profile.R` (`wock-schedbias-profile.rds`).

## 6. What is not yet known

1. **The bias at production settings.** Every bias measurement above uses a
   truth taken from a 96-compartment fit, and none was made at 384
   compartments or with the current priors together. A simulation study
   (SLURM 30875, started 15 October) measures it: two truths (the original,
   and the converged 384-compartment fit with per-trial cycle lengths
   39.4–45.6 h), simulated and fitted at 192 and 384 compartments crossed,
   true synchrony at 0.5, 1, and 3 times the fitted value, synchrony
   estimated or fixed at the truth, five datasets per cell. Reading rule,
   fixed before any fit: a cell needs three datasets passing the gate.
   **Partial read (96 of 240 fits, `simstudy-read.R`,
   `simstudy-read-2026-10-18-partial.txt`): no cell readable yet.**
2. **Fixing synchrony at 384 compartments does not sample.** In the same
   partial read, 23 of 24 fits with synchrony fixed and 384 compartments
   fail the gate, 20 severely, including fits with synchrony fixed at the
   true value on data the model generated itself
   (`simstudy-health-2026-10-18.txt`). So the fixed-synchrony diagnostic in
   section 2 cannot be repeated at 384 compartments, and the remaining fits
   of that arm were cancelled. The comparison of estimated against fixed
   synchrony will be readable only at 192.
3. **Some trials have two likelihood peaks.** With synchrony pinned at 15,
   the trial MMV048_PIB has profile-likelihood peaks at 38.5 h and 43.0 h,
   0.65 log-likelihood units apart (`cl-profile.R`,
   `cl-profile-read-2026-10-14.txt`). This is uncertainty rather than bias,
   and moving one trial by 4.5 h moves the 13-trial mean by at most 0.35 h.
4. **The real-data link between sampling density and cycle length.** In the
   96-compartment fits, trials whose series were sampled more sparsely had
   longer estimated cycles (correlation across 13 trials −0.93 on posterior
   means; −0.52 within draws, 95% interval −0.85 to +0.14). Simulating every
   trial from one cycle length reproduces a correlation of −0.19 to −0.89
   (7 datasets), so the design manufactures part of it, but none of the 7
   reaches the real value (`wockner-cl-clustering.R`,
   `schedsim-analyze.log`). Not re-examined at higher compartment counts.

## Consequences for reported numbers

- The current real-data estimates, mean over 13 trials at 384 compartments,
  are **41.34 h** with synchrony pinned at 15 (41.29–41.41 h over five runs,
  one of which passes the gate) and **41.82 h** with synchrony estimated
  (41.79 and 41.82 h in the two of four runs that pass)
  (`seed-probe-read.R`, `seed-probe-read-2026-10-15.txt`,
  `pin15-pooled-2026-10-14.txt`). Neither is corrected for bias, because the
  bias has not been measured at these settings.
- Any cycle length from this model is conditional on the compartment count
  and on how synchrony is handled, and should be reported with that
  sensitivity rather than as a single number.
- Correcting the real-data estimates with the 96- and 192-compartment biases
  places them in the low 40s (40.5–42.3 h; `nc-bias-correct.R`,
  `nc-bias-correct-2026-10-07.txt`). Those biases were measured under
  different priors and from a truth near 45 h, so this range indicates
  direction and is not an estimate.

## Sources

| claim | script (`_scripts/`) | saved output (`_data/`) |
|---|---|---|
| baseline simulated bias; intervention table | `wockner-schedule-sim.R`, `wockner-schedule-sim-analyze.R` | `schedsim-analyze-2026-10-05.txt`, `schedsim-analyze-2026-10-06.txt`, `schedsim-analyze.log` |
| likelihood share of the bias | `phase-volume-test.R` | `phase-volume-test.log` |
| per-trial ML bias near zero | `wockner-schedule-bias-profile.R` | `wock-schedbias-profile.rds` |
| summary statistic ruled out | `summary-stat-check.R`, `wockner-schedule-sim-mode.R` | `summary-stat-check.log`, `wock-schedsim-mode.rds` |
| bias under corrected priors; synchrony attenuation | `schedsim-truthfit-check.R` | `schedsim-truthfit-check.log` |
| real-data prior sensitivity | `wockner-prior-panel.R`, `panel-bshape-spread.R` | `prior-panel.log`, `panel-bshape-spread.log` |
| compartment-count ladder on real data | `nc-ladder-read.R`, `nc-ladder-loo.R`, `nc-numerical-check.R` | `nc-ladder-2026-10-06.txt`, `nc-ladder-loo-2026-10-18.txt`, `nc-numerical-check-2026-10-06.txt` |
| compartment count without priors | `decay-law-read.R` | `decay-law-read-2026-10-08.txt` |
| observable period vs parameter | `nc-period-check.R`, `nc-mechanism.R` | `nc-period-check-2026-10-07.txt`, `nc-mechanism-2026-10-07.txt` |
| compartment-count misspecification | `nc-2x2-read.R`, `nc-bias-correct.R` | `nc-2x2-2026-10-07.txt`, `nc-bias-correct-2026-10-07.txt` |
| total spread fixed; profile over compartment count | `nc-profile-fast.R`, `nc-profile-fast-read.R` | `nc-profile-fast-2026-10-09.txt` |
| simulation study, partial | `simstudy.R`, `simstudy-read.R`, `simstudy-health.R` | `simstudy-read-2026-10-18-partial.txt`, `simstudy-health-2026-10-18.txt` |
| two-peak trial | `cl-profile.R`, `cl-profile-read.R` | `cl-profile-read-2026-10-14.txt` |
| current real-data estimates | `seed-probe-read.R` | `seed-probe-read-2026-10-15.txt`, `pin15-pooled-2026-10-14.txt` |

Full working record: `claude/findings.md` in the same repository.
