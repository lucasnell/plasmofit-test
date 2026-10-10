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
  - `claude/ipm-decision.md` — the age-structure design decision: chain plus
    a parameter, or an IPM, with the floor argument that settles it

## How this project works

- **Cluster, SLURM.** R is at `/programs/R-4.6.1/bin/Rscript`; the user library
  is only on the login-shell path, so run `bash -lc 'Rscript ...'` or the full
  path. Jobs are array scripts in `_scripts/*.sh`, output to `_data/`.
- **Never take more than half the node.** `cbsugreischar` is a single shared
  node: **256 CPUs and 1,031,340 MB**, so the budget is **128 CPUs and 515,670
  MB**. It is shared with other people's work — at the time of writing 250 of
  the 256 CPUs and 800 G were allocated to someone else — and the rule is
  about our own footprint, not about what happens to be free.
  **If a job, or our jobs taken together, would exceed either half, throttle
  the array**: write `#SBATCH --array=1-N%M`, where `N` is the number of
  array tasks and `M` is how many may run at once:

  ```
  M = floor( min( 128 / cpus-per-task , 515670 / mem-per-task-in-MiB ) )
  ```

  `%M` caps **simultaneously running** tasks, not the total, so all `N` still
  run — just fewer at a time. `M` for the profiles this project uses:

  | `--cpus-per-task` | `--mem` | CPU cap | memory cap | **M** |
  |---|---|---|---|---|
  | 1 | 8G | 128 | 62 | **62** |
  | 1 | 16G | 128 | 31 | **31** |
  | 4 | 16G | 32 | 31 | **31** |
  | 4 | 24G | 32 | 20 | **20** |
  | 4 | 48G | 32 | 10 | **10** |

  No script in `_scripts/` currently needs a throttle — the largest array is
  `wockner-schedule-sim-batch.sh` at 16 tasks × 4 CPUs × 16G = 64 CPUs and
  262,144 MB, half the budget on memory and half on CPUs. **But the rule binds
  on the aggregate of everything of ours running at once**, so check `squeue
  -u lan68` before submitting a second array: that script plus a 14-task 16G
  array would come to 504,832 MB, inside the limit by 2%.
  Re-derive the node figures with `sinfo -N -o "%N %c %m"` rather than trusting
  these if anything looks off. **Disk is not a node resource here** — `/home2`
  is shared Lustre (755 T, 201 T free) and `--array` throttling does not
  govern it; `_data/` is 5.2 GB and is the thing to watch.
- **Never edit a script while a job is reading it.** This has already produced
  one silently invalid result (group factor levels reordered mid-run). Check
  `squeue` first.
- **Convergence gate is max R-hat < 1.05.** A run that fails it is not a
  posterior and is not reported, even when most chains agree. Check `lp__` per
  chain before diagnosing: one stuck chain and a general failure want
  different fixes.
- **Every table of numbers must say what is in its cells.** A number here is
  never self-describing: `8.89` could be a posterior mean, a median, one
  group's value or a mean over groups, and `−40%` could be relative to a
  truth, to another arm, or to a prior. State it once, immediately above the
  table, covering as applicable: **what the number is** (posterior mean,
  median, mean over groups, sum, difference, ratio); **its units** (hours,
  log10 units, percent of what); **what it is relative to**, and which
  direction is better or less biased; **what it is aggregated over** (groups,
  trials, replicates, draws) and how many; and **whether it is paired**, since
  paired and unpaired numbers of the same quantity differ here by more than
  the effects being measured. The column header is not enough — it names the
  quantity but not how it was computed, and the whole value of these notes is
  that a number can be re-derived a month later. The same applies to a number
  quoted in prose: write "posterior mean over 14 `grp_init` groups", not "the
  estimate", and name the script and saved output it came from. Purely
  descriptive tables — the script list, the `_data/` naming key — are exempt,
  since their cells are prose.
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

- **2026-10-08 — The √ decay law is not established, only unbeaten.** The
  decay-law test (SLURM 29684) could not distinguish synchrony decaying as
  √(cycles) from a component linear in cycles: at the pre-registered sigma =
  0.033 a linear component buys at most +0.396 and typically +0.0000
  log-likelihood units over 1130 observations. **An IPM therefore cannot be
  justified by appeal to the decay law**, and a 1-D Gaussian-kernel IPM would
  reproduce √ by construction. A by-product of the same job: model A is a
  priors-free likelihood profile over `n_c`, and 192 beats 96 in **14 of 14**
  units by **+72.4** summed log-likelihood units, so **`n_c` = 96 being wrong
  is not a prior artefact**.

- **2026-10-09 — The data constrain the TOTAL stage spread, not how it splits
  between initial synchrony and accumulated desynchronisation.** Across rungs
  where the accumulated spread varies 2.0x (chain) and 4.4x (gamma IPM), the
  total spread at the last observation varies by only 12.7% and 23.8%;
  `b_shape` falls monotonically as the rung rises in **14 of 14 units**
  (Spearman −1.00). This explains the whole age-structure thread: why the
  decay law was not learnable, why the `n_c` profile saturates, and why
  `b_shape` is unidentified. At the best rungs the split is ~0.11 cycles
  initial against 0.05–0.08 accumulated, so **the desynchronisation is mostly
  inherited from t = 0, not acquired during the window.**
- **2026-10-09 — Do not build the IPM.** Both profiles SATURATE: the data put
  an upper bound on dispersion and no lower bound, so an IPM would return
  `sigma_d` against zero — a boundary estimate, not a rate. The cost argument
  never applied: the convolution that was supposed to supply it is 5x SLOWER
  than the production series, and was removed from the package again on
  2026-10-09. Non-uniform stage rates do not help either,
  since they can only ADD dispersion above the Erlang floor and the data want
  less, not more.
- **2026-10-09 — `cycle_length` is stable once `n_c` >= 384**, spanning 40.93
  to 41.15 h across the 2-log-likelihood interval of the ML profile. The 3.20 h
  span of the original ladder came from `n_c` = 96 and 192, now known wrong.
  The first converged production fit (`n_c` = 192, `b_shape` 400, SLURM 30525)
  gives **41.20 h**, mean over 13 trials, with 5.8% divergences. **The second
  rung (`n_c` = 384) has failed twice** — 29635 and the reseed 30829, the
  second with all four chains apart (`lp__` spread 169.8, max R-hat 7.74) —
  so the ladder is not closed and no cycle length is final.

## Working notes

- Current status: `PROJECT_INDEX.md`
- Live work: `TODO.md`
- Last session: `handoff.md`
