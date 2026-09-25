# Handoff: the `b_shape` ladder and three refits, 2026-09-25

*Transient. Delete once the results are in `findings.md` and threads 2 and 3
are updated.*

## What is running

| job | what | tasks | outputs |
|---|---|---|---|
| **28963** | `b_shape` ladder, `wockner-fit.R` 14-18 | 5 | `_data/wock-fit-np_bs{50,84,100,150,250}.rds` + `-LOO-` |
| **28964** | convergence refits on a new seed | 3 (4, 5, 29) | `_data/wock-schedsim-{fit,RES}-<arm>-rep<n>-seed1618033989.rds` |

Descriptive table, no numeric cells. Both mail `lan68@cornell.edu` on
`END,FAIL`. Real-data fits are ~1.5-2 h; the simulation refits ~2 h.

## 28963: the `b_shape` ladder

`b_shape` is the shape of a symmetric Beta over cycle position, so it is a
**synchrony** parameter and converts to a starting-stage age range (see
`findings.md`, "What `b_shape` means biologically"). The default prior's 14.9
implies a **20.4 h** spread on a ~45 h cycle, which is nearly no synchrony.
`mmcm.pdf` Fig. 1 assumes ~9 h for controlled human infection trials, which
inverts to `b_shape` ~ 84.

Five rungs, `b_shape` pinned by a tight prior (`sd_log_b_shape = 0.05`, ±10%
at 95%) at **50, 84, 100, 150, 250**, all carrying the corrected
`log10_total0` prior (`sd_log10_total0 = 1`). So the ladder reads against
**task 9, `np_wide_total0`**, not against `no_pool`.

`max_shape` is raised 250 -> 400 for all five, so the 250 rung is not sitting
on its own boundary. Every rung's prior is far from 400. This is the one
respect in which tasks 14-18 differ from tasks 1-13.

**A ladder, not one fixed value, on purpose.** Fixing `b_shape` converts an
uncertainty into an assumption; the sensitivity across the ladder is the
output, not the value at any one rung.

### How to read it

1. **`cycle_length` across the rungs.** The prediction recorded before the
   runs: these should move it only **~0.1-0.3 h** relative to
   `np_wide_total0`, because 65 -> 100 narrows the age range 10.1 -> 8.2 h
   where 15 -> 65 narrowed it 20.4 -> 10.1 h. **A much larger move means that
   reasoning is wrong** and the ladder, not the prediction, is right.
2. **elpd across the rungs**, against `np_wide_total0` (free `b_shape`) and
   `np_wide_both`. Flat elpd means fixing anywhere in the ladder is
   defensible and the range is what gets reported. Elpd falling off sharply
   above ~84 means the data object to the high-synchrony assumption, which is
   worth knowing and should be reported rather than smoothed over.
3. **Sampler health.** The widened-`b_shape` fits carry the batch's worst
   divergences (3.4-4.0%). Pinning a bounded scalar should sample *better*,
   not worse. If a rung is worse, that is a finding about the geometry, not a
   nuisance -- and it is not the `tight_sigma` trap, which pinned a
   hierarchical scale near zero and made a funnel.

Read all three together. A rung that fits as well, samples better, and is
biologically defensible is the case for fixing `b_shape`; a rung that buys
biology at a real predictive cost is a tradeoff to state, not to resolve
silently.

## 28964: the three refits

`default-rep4` (max R-hat 1.22, 478 divergences, min ESS 13),
`default-rep5` (1.07), and `wide_total0-rep4` (1.06) failed the gate in
SLURM 28940. Per thread 9, both previous failures of this kind were bad fit
seeds rather than hard datasets, so this changes `SCHEDSIM_FIT_SEED` to
**1618033989** and nothing else. Output names carry the seed and the failed
runs stay on disk beside them.

**Why it matters more than two lost arms**: `default` is the baseline every
paired contrast in the cycle-length budget subtracts, so with reps 4-5
missing the previous batch added no *paired* replicates at all and the budget
is still n=3.

**If a refit fails again**, that is evidence the dataset is hard rather than
the seed, which would be new -- thread 9 found the opposite both times. Do
not simply try a third seed without recording that.

## What to do when they land

1. `squeue -u lan68`; then `sacct -j 28963,28964 --format=JobID,State,ExitCode -P`.
2. Ladder: read the fits directly, or extend `_scripts/wockner-prior-panel.R`'s
   `PANEL` to include the five rungs -- it already does sampler health,
   posterior means, `loo`, and paired contrasts, and a missing fit becomes
   `NA` rows rather than an error.
3. Refits: re-run `_scripts/wockner-schedule-sim-analyze.R`. It now **errors**
   if two converged `default` runs exist for one replicate, rather than
   silently duplicating rows through the join. If that fires, decide which
   seed is the baseline before pairing.
4. Write into `findings.md` with `Cells:` per `conventions.md`; update
   `CLAUDE.md` and `threads.md`; delete this file.

## Still not started

- **Option 4, thread 4**: `hold_out` masking then Design A. The only honest
  route to the hierarchy answer. All-zero `hold_out` reproducing current fits
  exactly is the regression test.
- **Option 6**: cycle-length bound geometry, `[35, 50]` *moved* not widened.
  Now the one cheap untested suspect for the ~1.5 h of cycle-length bias that
  survives the corrected priors.

## Carried over, unchanged

- The package default `sd_log10_total0 = 0.25` is still in place by decision,
  and `wockner-fit.R` still inherits the package defaults rather than pinning
  them -- untested because `np_anchor` passes `inoc_size`, which replaces
  that prior. Test before pinning.
- `claude/handoff.md` and `claude/handoff-2026-09-25.md` are both spent and
  awaiting a decision to delete; their content is in `findings.md`.
