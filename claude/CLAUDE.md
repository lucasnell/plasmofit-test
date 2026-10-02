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

## Resume here, 2026-10-02 (end of session)

**26 tasks running**, all mailing `lan68@cornell.edu` on `END,FAIL`. Nothing
below waits on anything else. Both repos are committed; a restart costs only
the wait.

| job | n | what |
|---|---|---|
| 29537 | 1 | thread 6 drift, `np_wide_total0-drift` (pre-change build) |
| 29545 | 1 | thread 6 drift, `np_bs400-drift` — makes it n=2 |
| 29538 | 4 | **Design A**, configs 23–26, last-third mask |
| 29543 | 16 | budget replicates 6–7 (6 arms), bound arms reps 4–5 |
| 29544 | 4 | third-seed refits: `default-rep4`, `cl_move-rep2`, `cl_wide_move-rep2/3` |

**Do not edit `_scripts/wockner-fit.R` or `_scripts/wockner-schedule-sim.R`
until the queue drains.** Both are being read; `Rscript` parses incrementally
and an edit mid-run has killed jobs here before.

### First three things to do, in order

1. **Submit Design A at a quarter-mask** (4 fits) the moment the queue
   drains. It needs 4 configs in `wockner-fit.R` mirroring 23–26 with the
   mask at a quarter instead of a third — the `ho_last_third` column is built
   in that script, so add an `ho_last_quarter` beside it. Running both
   fractions means the fraction cannot be chosen after seeing the answer,
   which is the whole point of having pre-registered it.
2. **Read the five running jobs back** against the decision rules below.
3. **Then the `sd_iRBC` arm** — see "where the bias question stands".

### Decision rules, fixed in advance

Only a result falling *outside* one of these is worth stopping to discuss.

1. **Thread 6 drift (29537, 29545).** Add both `-drift` pairs to `PAIRS` in
   `_scripts/wockner-bshape-regression.R` and re-run. Read mean |z| against
   the seed-only null's **0.85** and the b_shape test's **1.23**. Near 1.23 →
   the excess is rebuild drift, the b_shape change is exonerated, the merge
   stands, thread 6 becomes a reporting caveat on `R`. Near 0.85 → the excess
   belongs to the change; revisit the merge before anything else.
2. **Design A (29538).** Score **held-out elpd only** — `log_lik` summed over
   `hold_out == 1`. `generated quantities` computes it for every observation
   fitted or not, so the total mixes the two and answers nothing. The mask is
   in `_data/wock-data-daA_*.rds`. Read the quarter-mask run beside it before
   writing a verdict.
3. **Refits (29544).** If `default-rep4` fails a **third** seed, stop
   reseeding and record it — that contradicts thread 9 and is a result, not
   an obstacle.
4. **Budget replicates (29543).** Re-run
   `_scripts/wockner-schedule-sim-analyze.R`; it globs, no argument needed.
   Expect the `b_shape` prior's share to move again (0.50 → 0.43 h from n=3
   to n=4). Quote the n=6–7 figure and its range.

### Where the bias question stands

The ~1.5 h of unexplained cycle-length bias is **not fit-limited** — none of
the 26 running jobs touches it. This session closed two candidate remedies
and left one live lead:

- **The posterior median does not help** and moves the bias *up* by
  0.03–0.09 h; the per-trial posteriors are left-skewed. Mean, median, and
  mode span 0.05 h on a 1.8 h bias. It *does* matter for `b_shape` (~10%),
  whose median should be what gets quoted.
- **The MAP cannot settle it either.** The joint surface has modes spanning
  ~2000 nats and 7–13 h of `cycle_length`; the best was never found twice in
  ~200 starts and more searching keeps finding better ones.
- **Still live and never isolated**: the MLE fixes `sd_iRBC` at truth where
  the fit estimates it. That is the last named component of the
  MLE-to-posterior gap that nobody has separated, and it is **one simulation
  arm** — add `sd_iRBC` fixed at truth to `ARMS` in
  `wockner-schedule-sim.R`, paired against `default`. Cheapest untried thing
  on the question.
- Beyond that it is idea-limited, not compute-limited.

### Package repo state

`main` at **`5be1b47`**, **not pushed**. Carries the merged `b_shape`-as-data
change and `hold_out` masking on top. The installed build is this one.

Five files show as modified and are **deliberately uncommitted build noise**:
`R/RcppExports.R`, `src/RcppExports.cpp`, `src/available.stanfunctions.cpp`,
and `stanExports_archer_fit_single.h` / `stanExports_test_ode.h` — all
namespace-hash churn from the rebuild, with zero `hold_out` references. The
four model headers that *did* change in substance are committed. Do not
"tidy" the rest into a commit; they regenerate on every install.

### Scratch to clean up

`/home2/lan68/plasmofit/.prechange/` holds a git worktree at `e61fb47` and a
library built from it, used only by 29537 and 29545. Once thread 6 is written
up: `git worktree remove /home2/lan68/plasmofit/.prechange/src` and delete
the directory.

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
