# Handoff — 2026-10-09

## State of play right now — read this first

**One job running: SLURM 30829**, the entry 40 reseed, submitted 2026-10-09,
**~24.5 h**. It is the last rung of the production ladder and the only thing
blocking a reportable cycle length. Everything else is finished and read:
30525 (the production reseed), 30641 + 30823 (the eight-start dispersion
profile), 30524 (`n_c` = 768) and 30527 (`n_c` 128/256/512).

**The age-structure question is answered**, and the IPM is closed. See
"The job that is running" for how to read 30829 when it lands.

```bash
cd /home2/lan68/plasmofit/plasmofit-test
squeue -u lan68          # one job: 30829
```

## The answer: the data constrain the TOTAL spread, not its split

This is the result that closes the thread, and everything else follows from
it. Fitted stage-distribution spread in **cycles** at the last observation
(k = 4.8), pooled over the 14 `grp_init` units: `sd_init` from the fitted
`b_shape`, `sd_acc` = √(4.8/rung) from the rung, `sd_tot` their quadrature sum.

| profile | accumulated spread varies by | **total** spread varies by |
|---|---|---|
| exact chain, 9 rungs | **2.0×** | **12.7%** of its mean |
| gamma IPM, 16 rungs | **4.4×** | 23.8% of its mean |

**The fit holds the total nearly fixed and trades initial against
accumulated.** `b_shape` falls monotonically as the rung rises in **14 of 14
units** — Spearman correlation of log(`b_shape`) with log(rung) is **−1.00**.

That one fact explains four separate open questions at once:

- **why the decay law was not learnable** (29684) — if the end-of-window
  spread is nearly fixed, its *growth law* is barely constrained;
- **why the `n_c` profile saturates** — any rung fits, because `b_shape`
  compensates;
- **why an IPM's `sigma_d` would sit at a boundary** — the data constrain the
  sum, not the accumulation term;
- **why `b_shape` is unidentified** — observations begin 1.6 cycles in, so it
  and the accumulation rate enter only through their sum.

At the best rungs the split is ~0.11 cycles initial against 0.05–0.08
accumulated: **the desynchronisation is mostly inherited from t = 0, not
acquired during the observation window.**

## What follows: do not build the IPM

Both profiles **SATURATE** under the rule pre-registered in
`nc-dispersion-profile-read.R`: the data put an upper bound on dispersion and
**no lower bound**, so they do not exclude zero. An IPM would return `sigma_d`
pressed against zero with an interval touching it — a boundary estimate, not
a rate.

| profile | best | 2-unit interval | top-rung gain | unimodality violation | verdict |
|---|---|---|---|---|---|
| exact chain, 64–1024 | 1024 | {512, 1024} | +1.30 | 0.00 | SATURATES |
| gamma IPM, mesh fixed 192 | 2264 | {1683, 4096} | +0.69 | 0.69 | SATURATES |

Non-uniform stage rates do not rescue it either: they can only **add**
dispersion above the Erlang floor, and the data want less, not more. And the
cost argument had already gone, when the exact chain turned out to be a
convolution that is 38× faster with no rewrite and no change of model.

`claude/ipm-decision.md` holds the full case; this is its conclusion.

## The cycle length, and why it is not final

**SLURM 30525 converged.** Entry 39 (`np_bs400_nc192`), reseeded: **max R-hat
1.0432**, inside the gate, with `lp__` means −834.7, −828.3, −836.5, −830.1 —
a spread of **8.3** against 29635's ~96. The stuck chain is gone; reseeding
was the right diagnosis and the tuned retry was not needed.

**Mean `cycle_length` = 41.203 h** over the 13 trials, range 40.762–41.982, sd
across trials 0.412.

Two caveats that travel with it:

- **5.8% divergences** (233 of 4000), above the 4.0% that was previously the
  worst in the project. This is a mean over an imperfectly explored posterior.
- **It is one rung.** Entry 40 (`n_c` = 384) failed in 29635 and has never
  been reseeded, so the question the ladder exists to answer — does the `n_c`
  effect survive a pinned `b_shape`? — is still open.

Independently, the ML profile says `cycle_length` is **stable at 40.93–41.15 h
once `n_c` ≥ 384** (spread 0.22 h across the 2-log-likelihood interval),
against a 3.20 h span over the original ladder. The large span came from `n_c`
= 96 and 192, now known to be wrong.

## The job that is running: SLURM 30829, entry 40

`_scripts/wockner-fit-nc-bs400-reseed40.sh`. `np_bs400_nc384`,
`WOCKFIT_SEED=20261009`, `WOCKFIT_SUFFIX=-seed2`, walltime 3 days against a
measured 24:34:15.

**Why it matters.** Entry 39 converged, so one rung of the production ladder
exists and the other never has. The question the ladder exists to answer —
does the `n_c` effect on `cycle_length` survive a pinned `b_shape`? — needs
both and cannot be answered from one.

**How to read it.**

1. **`lp__` per chain first**, before the R-hat gate. One stuck chain and a
   general failure want different fixes, and the stuck chain is this config's
   known failure mode (29635_40: lp__ −996.0 against −925.6, −955.3, −957.8).
2. **Max R-hat < 1.05.** A run that fails is not a posterior and is not
   reported, however well the other chains agree.
3. **Compare `cycle_length` against entry 39's 41.203 h.** The deterministic
   tests predict the `n_c` effect **persists and grows** with `b_shape`
   pinned: 0.98 h at `b_shape` 400 against 0.72 h at 15, in all twelve cells
   of the least-squares table. **If it converges and the effect is gone, the
   mechanistic account is wrong** and should be revisited, not patched.
4. For context, the ML profile — `b_shape` free, no priors — puts
   `cycle_length` at **40.93–41.15 h once `n_c` ≥ 384**.

**If a chain sticks again** at a similar lp gap, the mode is real: report it
as multimodality, show both modes, and **do not reseed a third time**. The
tuned retry (`wockner-fit-nc-bs400-retry.sh`) remains unsubmitted and is now
*less* justified than when written — it raises `max_treedepth` on the strength
of 46–61% saturation from 29635, but 30525 saturated in only 2% of
transitions with a maximum of 10.

## Mistakes made this session, so they are not repeated

Four, all mine, all in the instruments rather than the science. They are in
`claude/gotchas.md` in full.

1. **A two-sided regression check failed six tasks for improving.** It was
   written for a two-start run, where values should match exactly; with eight
   starts a higher log-likelihood is an improvement, not drift. Ask which
   direction a difference can legitimately go, and re-ask when the thing being
   compared changes.
2. **That check ran before `saveRDS`**, so stopping threw away ~2.5 h of
   completed fitting per task. **Write results first, then judge them.**
3. **A loess residual was used to measure optimiser noise.** It measures
   *curvature*. The eight-start chain profile is strictly monotone — zero sign
   reversals — and still scored 2.51. Replaced with a **unimodality**
   violation, which scores 0.00 on that same profile.
4. **A mesh-convergence check was mistaken for independent confirmation.** It
   refit the same rung at twice the mesh and moved it 0.50, which looked like
   the drop was real; both meshes ran the same optimiser from the same starts.
   A resolution check cannot detect optimiser error.

The near-miss worth remembering: the two-start profile appeared to **turn
over**, which is the branch that would have justified building the IPM. It was
one unit's optimiser failure worth **+11.29** log-likelihood units.

## Context for the next session

**Arbitrating optimiser disagreements.** The two-start `mat_exp_series` path
says the profile TURNS OVER at `n_c` = 512; the eight-start convolution says
SATURATES. There is nothing to adjudicate: the forward maps agree to 1e−11,
every rung-wise difference is ≥ 0, and the two-start path is **4.10 units
short at `n_c` = 768**, which alone flips its 512 → 768 step from +0.61 to
−2.81. More starts can only raise a maximum, so the higher value wins.

**A flag, not a contradiction.** The free `b_shape` in the ML screen settles
near **9–15** at the best rungs, against the production pin of **400**, whose
initial spread is six times smaller. The screen has no priors and no
hierarchy, so the two are not directly comparable — but the gap is large and
has not been explained.

**What did not work, so it is not retried.** Sequestration-grid discretisation
as the cause of the `n_c` shift (wrong sign, eight times too small); the
pre-registered reading rule for `nc-mechanism.R` (the shift shrinks with
window length, not grows); and three cost extrapolations, every one low — the
worst by 9× across rungs of the *same* model, because the optimiser's
iteration count grows with the rung.

**Everything in `_data/` is gitignored except the reader outputs**, which are
negated in `_data/.gitignore` so a number quoted in the notes can be
re-derived from the repo alone. The fits themselves — 5.2 GB — exist only on
the cluster filesystem.

**Where things live.** Root `CLAUDE.md` (stable context and settled
decisions), `PROJECT_INDEX.md` (status, workstreams, decision log), `TODO.md`
(actionable layer), this file. Detail in `claude/`: `findings.md` (results,
every table's cells defined), `gotchas.md` (**read before running anything**),
`scripts.md`, `threads.md`, `references.md`, and `ipm-decision.md` (the
age-structure design decision, now closed).
