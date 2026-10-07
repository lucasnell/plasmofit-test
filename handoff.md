# Handoff — 2026-10-07

## State of play right now — read this first

**One job running: SLURM 29684, the decay-law test.** 14 tasks, one per
`grp_init` unit, submitted 2026-10-07 13:20, **~6 h** (8 h possible). Nothing
else is queued.

```bash
cd /home2/lan68/plasmofit/plasmofit-test
squeue -u lan68                                  # empty means it finished
ls _data/decay-law-unit*.rds | wc -l             # expect 14
/programs/R-4.6.1/bin/Rscript --vanilla _scripts/decay-law-read.R
```

If fewer than 14 files exist, check `_data/decay-law-<task>.err` and
`sacct -j 29684`. A missing unit is not fatal — the reader works on whatever
is present and prints how many it found — but say so when reporting.

**Do not edit `_scripts/decay-law-test.R` while 29684 is running.**

### What 29684 decides, and the rule for reading it

Does synchrony decay as **√(cycles)** (what the Erlang chain and a 1-D
Gaussian-kernel IPM both give) or **linearly in cycles** (between-parasite
heterogeneity in cycle duration)? Model A is the chain alone at `n_c` in
{96, 192, 384}; model B adds a lognormal mixture over `cycle_length`, 7
Gauss-Hermite nodes, at `n_c` in {192, 384}. B nests A at sigma = 0, so they
are compared at a **matched 6 parameters**, A spending one on choosing `n_c`
exactly as B spends one on sigma. `d_ll = ll_B - ll_A`, positive favours
linear.

- **B clearly ahead and consistent in sign across units** → the decay has a
  linear component, and **a Gaussian-kernel IPM is the wrong target**,
  because it would reproduce the √ law it was meant to replace. The fix is
  rate heterogeneity: a mixture, or a 2-D kernel over (age, rate).
- **A ahead, or a tie** → the √ law is adequate, and an IPM becomes a
  reasonable way to decouple the dispersion rate from the numerical mesh.
- **Mixed signs across units** → the design cannot tell, and the choice has
  to be made on biology rather than on these data.

A couple of units of total `d_ll` is noise at this parameter count. Expect
sigma around **0.033** if the linear law holds — that is the between-parasite
CV in cycle duration reproducing `n_c` = 192's spread of 0.158 cycles at
k = 4.8.

**Partial signal already visible**, not the headline and not to be reported
as a result: model A's own choice of `n_c` is coming out 192 > 96 in every
task so far, with 192 and 384 effectively tied in task 1 (−31.219 vs
−31.251). That is the same ordering and plateau the Bayesian 71.8 elpd result
gave, reached by maximum likelihood with no priors — independent
corroboration that `n_c` = 96 is genuinely worse.

### Nothing else is running, and one thing is waiting

`_scripts/wockner-fit-nc-bs400-retry.sh` is **written and deliberately not
submitted**. It retries the production `n_c` ladder (`b_shape` pinned at 400)
after SLURM 29635 failed to converge. **Submit a reseed-only version of entry
41 first** — see "Open follow-ups" below for why the script as written is
badly costed.

**Until that question is settled, no reported cycle-length number should
change.** That is the single most important standing constraint right now.

## Session topic

Reading back two jobs on the `n_c` question, then working out **why**
`cycle_length` moves with `n_c` at all. The session ended with the mechanism
explained, one job dead on convergence, and a new test running that decides
whether a structural rewrite is worth doing. Housekeeping: both repos
committed and pushed, 494 MB of scratch deleted, and these four files created
from `new-project.sh`.

## Key decisions

- **`n_c` is a biological assumption, not a numerical setting.** 96 predicts
  71.8 elpd worse than 192 with `b_shape` free. Ruled out an arithmetic cause
  directly: `max_rel_diff` is 8e−13 at `n_c` = 96 and *grows* with `n_c`, so
  if better arithmetic were the gain the higher rungs would have to be more
  accurate, and they are slightly less.
- **The mechanism is a one-sided observation window meeting a broadening age
  distribution.** Sequestration starts at a fixed 18.58 h of a ~45 h cycle, so
  ~40% of the cycle is visible; the stage distribution broadens at a rate
  `n_c` fixes and is right-skewed; the centroid of the *visible* subpopulation
  drifts, and that reads as a period change. At `n_c` = 96 the observable
  period runs ~0.75 h short of the `cycle_length` parameter.
- **Gate any IPM rewrite on the decay law.** A 1-D IPM with a Gaussian
  development kernel reproduces the same √t law the Erlang chain gives, so it
  would decouple the dispersion rate from the mesh without testing the form.
  SLURM 29684 tests √ against linear first.
- **Do not report anything from SLURM 29635.** Max R-hat 6.13 and 8.34 from
  one stuck chain. Three chains agreed, and reporting that subset would be
  choosing the chains that give a tidy answer.

## Open follow-ups

- [ ] Read SLURM 29684 with `_scripts/decay-law-read.R`. **~6 h from
      2026-10-07 13:20**, measured: model A takes 38 min per unit (61 s, 459
      s, 1777 s at `n_c` 96/192/384) and model B costs ~7x the same `n_c`
      because each likelihood runs seven quadrature nodes. B at 384 is the
      dominant term and is extrapolated, so 8 h is possible. Reading rules
      are in the script and in `claude/threads.md` thread 15.
- [ ] Resubmit the production `n_c` ladder — **reseed-only first**, entry 41
      (`n_c` = 192) alone, ~9 h. The written retry
      (`_scripts/wockner-fit-nc-bs400-retry.sh`) adds `adapt_delta` 0.95 and
      `max_treedepth` 12, which costs 2.5–5× because the sampler already
      saturated treedepth in 46–61% of transitions; entry 42 would likely
      exceed its own 3-day walltime.
- [ ] Decide whether to keep `claude/CLAUDE.md` — see below.

## Context for the next session

**What did not work, so it is not retried.**

- *Sequestration-grid discretisation* looked like the obvious cause of the
  `n_c` shift, and the source even documents an O(1/n_c) delay in
  sequestration onset. It is **wrong**: the duty cycle *rises* with `n_c`
  (0.39765 → 0.40242), implying `cycle_length` should rise by +0.284 h against
  an observed −2.17 h. Wrong sign, eight times too small.
- *The pre-registered reading rule for `_scripts/nc-mechanism.R`* — separate
  the channels by whether the shift is constant or growing with window
  length — **did not work**. The shift *shrinks* in every cell, because a
  longer window pins the period harder as well as accumulating more spread,
  and the second effect dominates. `_scripts/nc-period-check.R` is the clean
  instrument for the period channel.
- *A per-leapfrog cost probe* under-predicted the production fit by 2×,
  because the leapfrog count rises too. Cost is per-gradient × gradients
  per iteration; a fixed-iteration probe measures only the first.

**Corrections made to earlier claims in this session.** The forward map
accounts for **half to two thirds** of the real `n_c` shift, not the quarter
first estimated from the period channel alone. Both figures are in
`claude/findings.md`; the later one supersedes.

**A fact worth carrying.** Observations begin at 72 h — **1.6 cycles in** —
and run to 4.8 cycles. Nothing is observed in the first cycle and a half. That
is why `b_shape`, which describes t = 0, is unidentified: it is pure backward
extrapolation. A dispersion rate would govern change *within* the window and
should be better identified.

**Uncertain.** Whether the `n_c` effect survives a pinned `b_shape` is the
open question, and the deterministic tests *predict it should persist and be
larger* (0.98 h at `b_shape` 400 against 0.72 h at 15, in all twelve cells of
the least-squares table). If the retry converges and the effect is gone, the
whole mechanistic account is wrong and should be revisited rather than
patched.

**Notes restructuring.** `CLAUDE.md`, `PROJECT_INDEX.md` and `TODO.md` were
created at the repo root from `new-project.sh`. The detailed files in
`claude/` are unchanged and still hold everything: `findings.md` (results),
`gotchas.md` (read before running anything), `scripts.md`, `threads.md`
(long-form behind `TODO.md`) and `references.md`.
**`claude/CLAUDE.md` and `claude/conventions.md` were deleted** (approved
2026-10-07, recoverable from git history). `claude/CLAUDE.md` was superseded —
stable content to the root `CLAUDE.md`, status to `PROJECT_INDEX.md`, resume
block to this file — and partly stale: its "Package state" section still said
the package was on branch `fixed-b-shape` and unpushed, untrue since
`8dde0c1`. Before deleting, two items that existed **only** there were
migrated to `findings.md`: the bound-asymmetry result and the "Known thin
spots" caveats, now "Standing caveats on the evidence". All eight distinctive
figures from its 2026-09-25 sections were confirmed present in `findings.md`.
`conventions.md` was folded into the root `CLAUDE.md` under "How this project
works".
