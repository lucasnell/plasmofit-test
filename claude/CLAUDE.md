# plasmofit-test — scripts and fits

Driver scripts and saved fits for the `plasmofit` package; not itself a
package. **Working on the cluster**, primary directory
`/home2/lan68/plasmofit/plasmofit-test`, package source
`/home2/lan68/plasmofit/plasmofit`. Those are **two separate git
repositories** and both need committing separately; package-side notes live
in `plasmofit/claude/` (`chat3/CLAUDE.md` for model internals,
`chat5/CLAUDE.md` for the inoculum-anchored prior). These scripts used to
live in `plasmofit/_testing/`, so old commits still name that path.

This file is an index and the current state. **Read the file you need:**

| file | when to read it |
|---|---|
| `claude/scripts.md` | what each script does, the `CONFIGS`/arm mechanics, `_data/` naming, how to run on the cluster |
| `claude/gotchas.md` | **before running or editing anything** — the traps that have cost hours each |
| `claude/findings.md` | the modelling results: hierarchy, pooling offset, schedule-bias simulation, nuisance priors |
| `claude/threads.md` | open threads 1–10, in priority order |
| `claude/conventions.md` | how to write in these files; **every numeric table must define its cells** |

## Resume here, 2026-10-05

**6 tasks running**: SLURM **29612**, the cycle-length prior's **location**
(arms `cl_center_truth` and `cl_center_low`, tasks 71–73 and 78–80), ~2 h,
mailing on `END,FAIL`.

**Do not edit `_scripts/wockner-schedule-sim.R` while it runs.**

### How to read 29612

Re-run `_scripts/wockner-schedule-sim-analyze.R`; it globs and both arms are
now in `PAIRED_ARMS`. Read it as a **slope**, not three separate contrasts:
posterior cycle-length bias against prior-mean displacement from the truth,
across three arms —

| arm | `cl_prior_center` | prior mean | displacement |
|---|---|---|---|
| `default` | 48 | 47.44 h | **+2.43 h** |
| `cl_center_truth` | = truth | 44.61 h | **−0.41 h** |
| `cl_center_low` | 42 | 42.09 h | **−2.93 h** |

That slope is the prior's weight, measured. **If it is near 0.15**, the prior
explains ~0.36 h of the ~1.5 h and thread 2 keeps most of its gap. **If it is
much larger**, the prior explains more than this project has believed since
the schedule-bias simulation, and the 0.15 h on record was an artefact of
testing width instead of location.

### The other open lead, from today

The fits are **systematically phase-lagged and the lag grows along each
series** — −0.15 to −0.41 cycles at the last observation, growing by ~1.3×
from first to last, i.e. 8–18 h of mistiming on a 45 h cycle. That is new,
and it means the cycle-length bias is not an isolated parameter sitting off:
the fit does not track the oscillation's timing. Any explanation of thread 2
now has to account for the lag as well as the long period. Nothing is
currently testing this.

### Still open, in priority order

1. **Read 29612** as above.
2. **Reseed `daH2_no_pool`** — one fit, closes the horizon ladder; it is the
   only rung that showed a signal and the only one that failed the gate.
   `sbatch --array=33 --export=ALL,WOCKFIT_SEED=1123581321,WOCKFIT_SUFFIX=-seed2 _scripts/wockner-fit.sh`
3. **The phase lag** — no instrument proposed yet.
4. **More null runs** only if the one unexplained regression figure matters.
5. **Thread 5, Design B** — still argued against.

## Settled, 2026-10-04

- **The hierarchy question is NOT settled.** Design A is mask-dependent:
  −3.27 (z −4.5) at the last third, −0.04 (z −0.1) at the last quarter. 89%
  of the third-mask deficit sits on points the quarter also holds out, so it
  is a training-set effect, not a horizon one. An earlier claim that the
  hierarchy earns its keep is **withdrawn**.
- **The Design A mask no longer reorders rows** (`rank(time)`, not
  `arrange`). Mask identical at 306 points; the rebuilt data list matches the
  September group level order again. `daA_*` and `daQ_*` are the boundary —
  they carry the old ordering and are internally consistent.
- **`sd_iRBC` can be supplied as data** (package `94b6d99`), which is what
  makes the `fix_sd` arm possible. Diagnostic only: it exists for simulation,
  where the generating value is known.

## Settled, 2026-10-03

- **The hierarchy question is NOT settled, and Design A says why.** At the
  last-third mask, collapsing `cycle_length` costs 3.27 held-out log units
  (z −4.5); at the pre-registered last-quarter mask, −0.04 (z −0.1). 89% of
  the third-mask deficit sits on points the quarter mask also holds out, so
  it is not a horizon effect — the two differ in what the models were
  TRAINED on, not what they were scored on. **An earlier claim here that the
  hierarchy earns its keep is withdrawn.** Running both fractions is what
  caught it.
- **There is no rebuild drift.** Identical source rebuilt and refitted at the
  same seed gives mean |z| 0.79 and 0.76 against a seed-only null of 0.85.
  Thread 6 is closed and the caveat that `R` carries irreproducibility beyond
  its MCSE is **withdrawn**.
- **The posterior median does not help the cycle-length bias** and moves it
  *up* 0.03–0.09 h; mean, median and mode span 0.05 h on a 1.8 h bias. It
  does matter for `b_shape` (~10%), whose median should be quoted.
- **The MAP cannot settle it either**: modes span ~2000 nats and 7–13 h of
  `cycle_length`, the best was never found twice in ~200 starts, and more
  searching keeps finding better ones.
- **The budget is not settling.** The `b_shape` prior's share reads 0.50
  (n=3) → 0.43 (n=4) → **0.513 (n=6)**, range widening to 0.64 h. The
  replicate spread stays comparable to the effect being measured.
- **Bound asymmetry is ruled out** (+0.033 h, n=4). The `[30, 60]` arm's
  +0.471 h is the prior-width confound named at submission, not bound
  distance.

## Package state, 2026-10-02 — `b_shape` can now be supplied as data

**The installed `plasmofit` is from branch `fixed-b-shape`, not `main`, and
is not pushed.** `archer_stan_data(b_shape = )` takes a single value recycled
over `grp_init` groups, or one per group; `NULL` (the default) estimates it
as before. In Stan, `b_shape_free` is zero-sized when the switch is on and
`b_shape` is a transformed parameter assembled from whichever source is
active, so it stays in the output and every script reading it works unchanged.
Package commits `0bf12b0` and `e637d8d`; this repo's follow-through is
`c213cbb`.

**The Stan models were recompiled**, so per `gotchas.md` nothing fitted from
here on is bit-comparable with the fits already on disk. That is fine for
pooling across noise replicates and not fine for a bit-exact regression test.

**The end-to-end regression has NOT been run.** The package test suite passes
(20 new assertions, sampling tests included) and a fixed fit returns
`b_shape` = 400.0 exactly, but nobody has yet checked that `b_shape = NULL`
reproduces the existing `np_wide_total0` posterior, or that `b_shape = 400`
reproduces `np_bs400`. Both are ~2 h fits and both are the real test.

## Nothing in flight as of 2026-10-02
## Nothing in flight as of 2026-10-02

The queue is empty. SLURM 29493 (ladder extension) and 29494 (bound geometry)
completed and are written up in `findings.md`.

**Do not edit a script while a job is reading it.** `Rscript` reads source
incrementally; an edit mid-run killed two 1.5 h jobs.

## Settled, 2026-10-02 — the ladder plateaus, and bound geometry is out

- **The `b_shape` ladder plateaus at ~400**, a starting-stage age range of
  about **4 h**. Still climbing 250 → 400 (+0.96 elpd, z 3.3), flat 400 → 600
  (−0.11, z −0.3). Pushing past it costs sampling (5.92% divergences at 600).
  `max_shape` on its own is worth +0.00 elpd, so the extended rungs are
  comparable with the earlier ones.
- **`cycle_length` is flat from 250 upward**, 43.6–43.7 h, against 44.0 h at
  `b_shape` 50 and 45.5 h with `b_shape` free. Anywhere above 250 gives the
  same answer to within 0.1 h.
- The data want **tighter** synchrony than the 9 h assumed for controlled
  human infection trials in `mmcm.pdf` (`b_shape` ~ 84); that rung is
  detectably worse than 100, 250, and 400.
- **Bound geometry does not carry the cycle-length bias.** Moving the window
  changes the paired bias by **+0.059 h** (centred, same width) and +0.03 h
  (moved and widened), against the ~−1.5 h that would be needed, with the
  sign wrong. **Thread 2 is out of cheap suspects for the remaining ~1.5 h.**
- Two rungs fail the R-hat < 1.05 gate and are void: `np_bs150` (1.0513) and
  the `max_shape` control (1.0550). Nothing rests on either.

## Settled, 2026-10-01 — fixing `b_shape` is a large predictive gain

- **Fixing `b_shape` beats leaving it free by +28 to +36 elpd**, at every
  rung from 50 to 250, and beats merely widening its prior (+24.4). 10–12
  standard errors. The case for fixing it does not rest on the biology.
- **Where it is fixed matters much less, but the data are not indifferent**:
  3.6 elpd across 84–250, with higher consistently and detectably better
  (250 beats 100 by 2.26, se 0.50). **The ladder is still climbing at its top
  rung**, so 250 is the top of what was run, not an optimum.
- **`cycle_length` is conditional on the pin**: 44.0 h at `b_shape` 50 down
  to 43.6 h at 250, and 1.5–1.9 h below the free-`b_shape` fit. That spread
  is the sensitivity to report. `R`, `log10_total0`, and `sd_iRBC` are flat
  across the ladder — only `cycle_length` tracks it.
- **Pinning samples worse**: divergences rise 2.05% → 4.58% with the pin and
  `np_bs150` sits at the R-hat gate. The best-fitting rung is the
  worst-sampling one.
- **The `b_shape` prior's share of the simulated cycle-length bias is 0.43 h,
  not 0.50 h**, and its spread is four times wider than three replicates
  suggested. The old figure's precision was luck.

## Settled, 2026-09-25 — `b_shape`, and what the trend misfit actually was

- **The within-series trend misfit is the `log10_total0` prior, not
  `b_shape`.** Widening `b_shape` moves the median within-series range
  1.66 → 1.73 against an observed 2.563 and leaves the posterior predictive
  p-value at 0; widening `log10_total0` moves it to 2.29 and `ppp` to 0.155.
  The hypothesis that `b_shape` ~ 65 was the model absorbing lack of fit is
  **tested and rejected**.
- **`b_shape` is not identified: this design cannot separate 15 from 65.**
  Under identical widened priors a truth of 14.9 returns 19.5–27.9 and a
  truth of 65.0 returns 32.2–37.4 — a factor of 4.36 in the truth becomes
  1.53 in the estimate. The real data report 65.4, *above* what the design
  returns when the truth really is 65, so the attenuation runs the opposite
  way from the artefact explanation. **No point estimate of `b_shape` is
  defensible**, and the calibration must not be inverted.
- **Widening the nuisance priors does not remove the cycle-length bias.**
  Simulating from the widened-prior truth and fitting with widened priors
  still leaves +1.83 h. A shift in `cycle_length` when a prior is widened is
  not bias removal — in simulation or on real data.
- Replicates 4–5 left the recovery numbers essentially unmoved (`wide_bshape`
  `b_shape` 23.9 → 24.3), so that story is not an artefact of three noise
  draws.

## Settled, 2026-09-25 — the eight-fit prior panel

`_scripts/wockner-prior-panel.R`, output `_data/prior-panel.log`. All eight
fits clear R-hat < 1.05; divergences run 1.05–4.0% and are the standing
caveat, worst in the two fits with both priors widened.

- **Thread 2 answered: the `b_shape` prior carries cycle length on real data
  too.** Widening it moves `cycle_length` **−0.47 h** in isolation (against
  −0.501 h in simulation) and **−1.29 h** on top of a corrected
  `log10_total0` prior, and gains **+24.8 elpd (se 2.2)**. The two nuisance
  priors are additive in elpd but **super-additive on `cycle_length`**:
  jointly −1.1 h against −0.27 h from the one-at-a-time contrasts. On real
  data this is a **shift, not a bias** — it bounds the priors' share of the
  reported 45.3 h from below and says nothing about the remainder.
- **`b_shape` is not identified in location.** Under `lognormal(2, 0.5)` its
  posterior is barely narrower than its prior (CV 0.433 against 0.533), so
  **the 14.9 reported throughout this project is a prior artefact**; widened,
  it reads 65.4 with a posterior sd of 50. Same standing as `R`.
- **Thread 3 answered: the hierarchy conclusion survives, its number does
  not.** `pooled_cl` fails to beat `no_pool` at all three prior settings
  (−1.38, −1.67, −1.68 elpd, se ~1), so correcting the misspecified prior
  does not rescue the hierarchy. But the pooling offset runs **−0.304,
  −0.110, −0.509 h** across those settings — a factor of 4.6, and
  non-monotone. **Quote −0.509 h**, the best-fitting setting; −0.304 h was
  drawn under a prior costing 104.7 elpd.

## Settled, 2026-09-24

Tables and derivations for all of these are in `claude/findings.md`.

- **The `normal(1, 0.25)` prior on `log10_total0` is misspecified.** Relaxing
  it gains **104.7 elpd (se 12.8)** and moves `R` from 6.2 into the 15–18
  range burst size implies. The inoculum anchor adds nothing over simply
  widening it (2.13 ± 2.91). **`R` is not identified by these data.**
- **The `b_shape` prior owns ~0.50 h of the simulated cycle-length bias and
  the `log10_total0` prior owns none**, so the `log10_total0`/`R` ridge is
  not the mechanism. ~1.07 h of ~1.97 h is still unexplained.
- Dead explanations: the simulated designs are *not* less informative (0.97
  of the real design's information about the oscillation), and truth rebuilt
  from posterior **draws** rather than the mean vector rescues nothing.
- The anchor switched off reproduces the pre-change posterior (width
  inflation 1.01×); both non-converged replicates were bad **fit seeds**.

## Known thin spots

- The real-data prior comparisons are n=1 by construction — one dataset, one
  fit per setting — so shifts get a Monte Carlo z, not a confidence interval.
  The panel's `mean_z` says a shift beats MCMC noise and nothing more.
- Divergences are non-zero in all eight panel fits and reach 4.0% in
  `pl_wide_both`, so every posterior mean in the panel is a mean over an
  imperfectly explored posterior.
- The pooling offset (−0.509 h) is ~0.4 of the within-fit posterior sd of one
  trial's `cycle_length` (1.28 h under `no_pool` with both priors corrected).
- Three noise realizations; three posterior-draw replicates, usable only
  after a refit on a second seed. **No cycle-length number is yet
  biological** — the panel bounds the priors' share from below and stops
  there.
