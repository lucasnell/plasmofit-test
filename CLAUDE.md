# plasmofit — malaria within-host dynamics — project context

## Project

- **Code:** plasmofit
- **Full name:** Within-host dynamics of *Plasmodium falciparum* from volunteer infection studies
- **My role:** analyst / lead on the modelling
- **Question:** What is the asexual blood-stage cycle length, and how much of
  the estimate is driven by modelling assumptions rather than by the data?

## People

| Who | Role | What they own |
|---|---|---|
| Lucas Nell (`lucasnell`, lucas@lucasnell.com) | analyst | the package, the fits, these notes |

<!-- Fill in collaborators. Greischar & Childs is cited as prior art in
claude/references.md, not as authorship. -->

## Where things live

- **Data:** `_data/wockner-cleaned.csv` — derived, already cleaned. 1130
  observations, 177 series, 13 trials. Everything else in `_data/` is output
  and is **gitignored**.
- **Code, two separate git repos, both needing separate commits:**
  - `/home2/lan68/plasmofit/plasmofit` — the R/Stan package
    (`github.com/lucasnell/plasmofit`)
  - `/home2/lan68/plasmofit/plasmofit-test` — this repo, driver scripts and
    saved fits (`github.com/lucasnell/plasmofit-test`)
- **Detailed notes**, which these four files index rather than replace:
  - `claude/findings.md` — the modelling results, with every table's cells defined
  - `claude/gotchas.md` — **read before running or editing anything**
  - `claude/scripts.md` — what each script does, arm/config mechanics, `_data/` naming
  - `claude/threads.md` — the long-form record behind `TODO.md`
  - `claude/references.md` — external papers that bear on the model
  - `claude/conventions.md` — how to write in these files

## How this project works

- **Cluster, SLURM.** R is at `/programs/R-4.6.1/bin/Rscript`; the user library
  is only on the login-shell path, so run `bash -lc 'Rscript ...'` or the full
  path. Jobs are array scripts in `_scripts/*.sh`, output to `_data/`.
- **Never edit a script while a job is reading it.** This has already produced
  one silently invalid result (group factor levels reordered mid-run). Check
  `squeue` first.
- **Convergence gate is max R-hat < 1.05.** A run that fails it is not a
  posterior and is not reported, even when most chains agree. Check `lp__` per
  chain before diagnosing: one stuck chain and a general failure want
  different fixes.
- **Every table of numbers must say what is in its cells** — see
  `claude/conventions.md`. A number here is never self-describing.
- **Grouping structure**, which governs almost every design decision:
  `grp_init` = trial × inoculum size (14 units), `grp_R` = `grp_cl` = trial
  (13), `grp_sd` = trial × cohort (27). **Within one `grp_init` unit every
  series shares every parameter**, so a unit is one trajectory with its
  observations as replicates.
- **Observations start at 72 h (1.6 cycles) and run to 216 h (4.8 cycles)**,
  4–8 per series. Nothing is observed in the first one and a half cycles,
  which is why any t = 0 quantity is backward extrapolation.
- **Masks use `rank(time)`, never `arrange()`.** Reordering rows changes the
  group factor levels `archer_stan_data()` derives and silently invalidates
  comparison against every earlier fit.
- Package changes need `rm src/*.o src/*.so` and `--preclean`, then
  **verify the installed binary**, not the source.

## Decisions that are settled

- **2026-09-24 — The `normal(1, 0.25)` prior on `log10_total0` is
  misspecified; it is relaxed to sd 1.** Relaxing gains 104.7 elpd (se 12.8)
  and moves `R` into the range burst size implies. Corollary: **`R` is not
  identified by these data**, and the inoculum anchor adds nothing over simply
  widening the prior.
- **2026-10-01 — `b_shape` is fixed rather than estimated.** Fixing beats
  estimating by +28 to +36 elpd at every rung, 10–12 standard errors. The case
  does not rest on the biology.
- **2026-10-02 — It is pinned at 400.** The ladder plateaus there and
  `cycle_length` is flat from 250 upward. **Any cycle-length number is
  conditional on this pin and must be reported with the ladder.**
- **2026-10-04 — The hierarchy question has no robust evidence either way.**
  Design A is mask-dependent (−3.27 at the last third, −0.04 at the last
  quarter) and the horizon ladder is flat once a non-converged fit is
  reseeded. An earlier claim that the hierarchy earns its keep is withdrawn.
- **2026-10-06 — `b_shape` is the nuisance whose estimation carries the
  simulated cycle-length bias.** Pinning it at the truth removes −0.92 h of
  +1.97 h; widening its prior removes −0.51 h; `sd_iRBC` and `log10_total0`
  remove none. Of the five nuisances only `b_shape` moves `cycle_length`.
- **2026-10-06 — `n_c` is a biological assumption, not just a numerical
  setting, and 96 is wrong.** It predicts 71.8 elpd worse than 192. The effect
  is structural, not arithmetic: `max_rel_diff` is 8e−13 at 96 and *grows*
  with `n_c`, so better arithmetic cannot be the cause.
- **2026-10-07 — The observable period is not the `cycle_length` parameter.**
  Sequestration begins at a fixed 18.58 h of a ~45 h cycle, so only ~40% of
  the cycle is visible; convolving that one-sided window with a broadening,
  right-skewed age distribution shifts the apparent period. At `n_c` = 96 the
  observable period runs ~0.75 h short of the parameter. **`n_c` sets both the
  numerical mesh and the desynchronisation rate, and the model has no separate
  parameter for the rate.**

## Working notes

- Current status: `PROJECT_INDEX.md`
- Live work: `TODO.md`
- Last session: `handoff.md`
