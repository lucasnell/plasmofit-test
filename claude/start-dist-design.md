# Design: a shared, one-sided starting age distribution

*Drafted 2026-11-09 with Lucas, after `findings.md`, "The inoculum's ages".
**Built, validated, installed LIVE, and first fitted (entry 55, SLURM
35287) on 2026-11-18.** Status and the decisions taken are at the end.*

## What is wrong with the current start

`beta_starts(n_c, b_shape, b_offset, total0)` puts a symmetric Beta(b, b)
over the WHOLE cycle and shifts it by `b_offset`, wrapping at the cycle end.
Two consequences the biology does not allow:

- **It straddles the cycle end.** With the mode at 1-5 h the left tail wraps
  into late schizonts. In `np_wide_total0_nc384_t` that puts 7-42% of the
  starting mass (median 20%) into sequestered compartments, up to 39% past
  30 h in DSM265. An injected inoculum of circulating, Glycerolyte-frozen
  iRBCs cannot contain those.
- **It is fitted 14 times.** Every unit got the same 3D7-V2 material, so
  `b_shape` and `b_offset` should not vary by trial and inoculum size. 28
  free parameters stand in for 2.

## What is pinned, and what is not

| | status | source |
|---|---|---|
| One distribution for all 14 units | **firm** | one donor unit, one bank, diluted (Cheng 1997; Woolley 2021) |
| No mass past ~23-27 h | **firm** | venous blood holds only circulating iRBCs; model's own curve leaves 1% of the age-0 share by 27.5 h |
| Rings only, so tighter still | **likely** | Glycerolyte freezing (ATCC 50110); inocula dosed on ring counts (Woolley) — but 3D7-V2's own ring % is unpublished |
| Width at most sd ~5.9 h | **firm as a ceiling** | ages of a fully asynchronous donor's circulating iRBCs |
| Width at least k_donor / `n_c` (cycles²) | **soft, and joint with `n_c`** | donor ran 3.5-5 cycles before collection |
| **The mode** | **NOT pinned** | see below |

**The mode is the open one.** Two leads, both weak:

- *Fever timing.* Blood was taken 6 h after onset of high fever, and fever
  follows mass schizont rupture (Kwiatkowski & Nowak 1991), which would put
  the rings near 6 h old. But falciparum fever is "usually irregular,
  showing no distinct periodicity", and these were primary infections in
  non-immune adults, so rupture timing cannot be read off fever onset.
- *Cheng et al.'s 39% dividing in the first 24 h.* A two-brood device fitted
  to five growth curves under a 48 h cycle, not a measurement
  (`findings.md`).

So the design should **constrain the support and share the parameters, and
leave the mode free** inside the window.

## Proposed parameterisation

Replace the wrapped symmetric Beta with a Beta on a sub-interval. Writing
age in cycles, a = t / `cycle_length`:

> start density on a in [0, w], Beta(alpha, beta) rescaled to that interval,
> with w = `start_window` / `cycle_length` a DATA constant.

Parameterised by mean and concentration rather than (alpha, beta), because
the two decisions are "where" and "how tight":

> mu in (0, 1) — mean age as a fraction of the window
> kappa > 2 — concentration; alpha = mu kappa, beta = (1 - mu) kappa

- `start_window` = 24 h as the default, a data constant so it can be moved.
- The old behaviour is `start_window` = `cycle_length`, `mu` = 0.5 + offset,
  `kappa` = 2 `b_shape`, so **the current model is a special case** and the
  switch can default to it.
- `mu` and `kappa` are SHARED across `grp_init`, where `b_shape` and
  `b_offset` are per-unit. That is the change that removes 26 parameters.
- No wrapping, so no mass in schizonts by construction, and the
  `generate_starts` split into circulating and sequestered compartments at
  t = 0 can be dropped for this variant (an injected inoculum is all
  circulating).

Scale, for choosing priors: a uniform start over [0, 24 h] has sd 6.9 h;
the circulating-curve limit is sd 5.9 h; the current fit's units sit at
3.3-6.1 h.

## How to add it without breaking anything

An integer data flag, `start_dist` (0 = current wrapped Beta, 1 = windowed
Beta), in all four Stan programs, plus `start_window`. With `start_dist` = 0
every existing config reproduces bit-for-bit, so the eight-fit panel and
every saved comparison stand. `archer_stan_data()` gains `start_dist` and
`start_window` with defaults 0 and `max_cl`.

Work, in order:
1. `plasmofit.stan`: `windowed_beta_starts()` beside `beta_starts()`; a
   unit test that `start_dist` = 0 reproduces the old vector exactly and
   that the two agree when the window is the whole cycle and mu = 0.5.
2. The four `archer_fit*.stan` programs: declare `mu`/`kappa` (size 0 when
   `start_dist` = 0), branch in `generate_starts`, priors.
3. `archer-fit.R`: arguments, validation, the data list.
4. Reinstall with `rm src/*.o src/*.so` and `--preclean`, then **verify the
   installed binary**, not the source (`gotchas.md`): check `fit@model_pars`
   and that a fit honours `start_dist`.
5. Regression: refit one existing config with `start_dist` = 0 and compare
   posteriors against the saved fit, scaled by Monte Carlo error, against a
   null of two seeds of identical code (`findings.md`, "Regression test for
   the anchor").

## Decisions for Lucas

1. **Window width.** 24 h (round, just past the 5% circulating point at
   24.3 h), 27 h (1% point), or `cycle_length` / 2?
2. **Shared or hierarchical.** One `mu`/`kappa` for all 14 units, or a
   hierarchy with a small between-unit sd to absorb dilution and the <= 60
   min between thaw and injection?
3. **Priors on `mu` and `kappa`.** Flat within the window, or something
   centred on the middle of the ring window?
4. **Does this go in the package, or into a scratch copy first?** The
   package is byte-identical to `8dde0c1` and nothing needs reinstalling
   right now; a reinstall recompiles the Stan models and invalidates the
   anchor regression null.


## Status, 2026-11-18: built in scratch and validated

**Decisions (Lucas, 2026-11-09 to 11-18).** Window **13 h**, from the 3D7
ring stage (0.30-0.31 of the cycle; `findings.md`, "How long the ring stage
is"), not the 24 h circulating bound first proposed. **Shared** across the
14 units. **Flat priors** on the two Beta shapes, `start_alpha`,
`start_beta` ~ Uniform(1, `start_shape_cap`), with the **cap = 3**; window and
cap are both data constants. A flat prior on the shapes is not flat on the
width: on a 13 h window a cap of 3 gives prior sd of the starting ages 2.87
h on average and P(sd > 3 h) = 0.27, against 3.75 h for a uniform start. Run
only the 13 h arm. Build in a scratch copy first.

**Where it lives.** Scratch clone `/home2/lan68/plasmofit/plasmofit-dev`,
branch `start-dist`, commit d39b20e on top of the live package's bba580f;
installed into the scratch library `/home2/lan68/plasmofit/.Rlib-dev` (every
other package symlinked from the live library). **The live package repo and
the live library are untouched.** Source-only diff:
`claude/start-dist-source.diff` (the commit also carries the stanExports
headers rstantools regenerates). Install log
`_data/install-start-dist-2026-11-18.log`.

**What changed.** New Stan functions `windowed_beta_starts`, `make_DC_win`,
`generate_starts_win`, `ew_poly_series_win`, `traj_combo_partial_sum_win`; the
four `archer_fit*` programs gain data `start_dist`, `start_window`,
`start_shape_cap` and parameters `start_alpha`, `start_beta`, and branch on
the flag. `archer_fit_single.stan` is unchanged. `archer_stan_data()` gains
the three arguments and checks them; `archer_fit()` fills them in for data
lists saved before they existed. **The old code path is untouched text**, a
separate partial-sum function rather than a branch inside the old one.

Two departures from the plan above:
- **5-point Gauss-Legendre**, not the 3-point of `beta_starts()`, so the
  shapes can go down to 1 (uniform) without the numerical guard setting the
  biological ceiling. Its worst error is 4e-5 in absolute compartment weight,
  for shapes strictly between 1 and 2.
- **The whole inoculum is circulating at t = 0** (added the same day at
  Lucas's request, commit 037fb6b). The old start splits each starting stage
  by y(age), the composition of a population that aged into it; an
  injected inoculum was drawn from venous blood, so it is all circulating
  and only later ageing sequesters it. In ew_poly_series()'s derivation that
  is N_j(0) = W_j, so D_d uses W_j / c_j and C_d (the total) is unchanged.
  **It has no effect on the fit to these data**: it raises predicted
  circulating parasites by ~2.7% for the first half-cycle and is identical
  to 8 decimals from ~36 h, because by then every injected parasite has
  burst once and its progeny enter circulation regardless. Observations
  start at 72 h, where the first cycle carries 3.5e-29 of the mass.

**Validation** (`_scripts/start-dist-validate.R`,
`_data/start-dist-validate-2026-11-18.txt`; no sampling, the compiled model
evaluated at fixed points):
1. `windowed_beta_starts()` against exact `pbeta` differences: machine
   precision at shapes 1, 2, 3; worst 3.9e-5 over n_c {96, 384, 512} x three
   cycle lengths x six shapes; continuous as the window edge crosses a
   compartment boundary.
2. **`start_dist` = 0 is identical to the old build**: on the production
   data (`np_wide_total0_nc384_t`), log_prob and its gradient at 8 fixed
   points differ from the live build by exactly 0 (904 gradient entries).
3. `start_dist` = 1, 13 h, cap 3, production data: 73 unconstrained
   parameters (113 minus 14 `b_shape_free` and 14 x 2 `b_off_vec`, plus 2);
   `b_shape` and `b_offset` zero-sized; log_prob and gradient finite; the
   gradient in `start_alpha`, `start_beta` matches central differences to
   8.8e-10.
4. The matrix_exp cross-check at n_c 96: max_rel_diff 6.0e-13 with the
   windowed start against 5.9e-13 for the old one, so the fast likelihood
   path and an independent matrix exponential of the same windowed starting
   state agree.
5. `wockner-fit.R` entry 55 builds the same data list as entry 54 apart from
   the three new fields, a missing window is refused, and the script's own
   `draw_as_list()` + `unconstrain_pars()` round trip works for both flags.

**Installed live 2026-11-18** from the scratch clone at commit **037fb6b**
(`R CMD INSTALL --preclean` into the live library
`/home/lan68/R/x86_64-pc-linux-gnu-library/4.6`; log
`_data/install-start-dist-live-2026-11-18.log`), and the INSTALLED binary
re-validated (`_data/start-dist-validate-live-2026-11-18.txt`): start_dist =
0 matches the old build's saved log density and gradient exactly, and the
windowed cross-check is 6.0e-13 with the all-circulating start. **The live
library therefore runs plasmofit-dev 037fb6b, which is NOT on the package
repo's `main` (still bba580f).** Porting the branch is Lucas's call. Do not
delete `/home2/lan68/plasmofit/plasmofit-dev`.

**First fit: entry 55, `np_start13_nc384_t`, SLURM 35287**, submitted
2026-11-18 (pinned, 40 CPUs; ~19 ms per gradient at the start).


## start_dist 2 and 3, 2026-11-19: testing why entry 55 failed

Entry 55 (`start_dist` 1) failed the gate with trials' `cycle_length` at
the bounds (`findings.md`, "Entry 55"). Hypothesis: sharing the start took
away the per-unit PHASE that entry 54's `b_offset` supplied. Lucas chose to
test it two ways ("do 1-2"):

- **`start_dist` 2: a mode per unit, one concentration shared.** A shift
  of a 13 h-wide shape inside a 13 h window is impossible, so "shared shape,
  per-unit shift" became per-unit MODE with shared CONCENTRATION.
- **`start_dist` 3: mode and concentration per unit** -- the window kept,
  the sharing dropped. With 2, this separates "shared" from "windowed".

**Parameterisation.** a = 1 + m (k - 2), b = 1 + (1 - m)(k - 2), with mode
m = `start_mode` in [0, 1] (a fraction of the window) and concentration
k = `start_conc` = a + b in [2, 2 `start_shape_cap`]: both shapes stay >= 1
for any m, the mode of Beta(a, b) is exactly m, and k = 2 is uniform. Flat
priors on m and k. The upper bound 2 x cap makes a symmetric start (m = 0.5)
carry the same cap on a and b as `start_dist` 1, but this is a different
prior from 1's flat prior on (a, b), not the same model with a constraint
lifted. Entry 55 put a at its floor of 1, i.e. m = 0; here m = 0 is inside
the allowed range rather than on a boundary of a shape parameter.

**Build.** Scratch clone commit **fc7db77** (branch `start-dist`), installed
ONLY into the scratch library; the live library and the package repo's
`main` stay at 13db2bd. `traj_combo_partial_sum_win` now takes per-group
shape vectors. Validation (`_scripts/start-dist-validate.R`, section 6;
`_data/start-dist-validate-23-2026-11-19.txt`):
- `start_dist` 0 identical to the old build (difference 0);
- `start_dist` 1 against the 13db2bd build at 8 fixed points: log_prob
  within 7.3e-12 (on values ~5e4) and gradient within 4.8e-16 -- rounding
  from the vectorised call, same model;
- 2: 86 parameters, 3: 99; the (m, k) -> (a, b) map exact (1e-16); finite
  log_prob and gradient; spot finite differences within 5e-7; the
  matrix_exp cross-check 6.6e-13 (2) and 1.0e-12 (3);
- `wockner-fit.R` entries 56 and 57 build entry 54's data plus the start
  fields, and the script's helpers round-trip both.

*Correction, 2026-11-20:* the script's `withS()` set a 1 h window for any
`start_dist` other than 1, so the log_prob, finite-difference and
matrix_exp numbers above for 2 and 3 were at a 1 h window, not 13 h. The
fits were unaffected: their data come from `archer_stan_data()` with 13 h,
which the last bullet checks. Rerun at 13 h
(`_data/start-dist-validate-below-delay-2026-11-20.txt`, 6b-6c): finite
differences within 3.9e-8, matrix_exp cross-check 6.5e-13 (2) and 9.9e-13
(3).

**Runs.** Entries 56 (`np_start13mode_nc384_t`, `start_dist` 2) and 57
(`np_start13unit_nc384_t`, `start_dist` 3), SLURM 35288, submitted
2026-11-19 with `WOCKFIT_LIB=/home2/lan68/plasmofit/.Rlib-dev`.


## Below age 0 and a delay, 2026-11-20: entries 58 and 59

**Why.** Entries 55-57 all fail the gate the same way, and 57's best region
fits 17-18 summed log-likelihood units below entry 54 with as many start
parameters. Entry 54's start mass outside [0, 13] h sits in the cycle's last
6 h (`findings.md`, "Entries 56 and 57"), so the data seem to want parasites
"younger than age 0". Two ways to give the start that phase (Lucas: "do
both"):

| | 58 `np_start13below6unit_nc384_t` | 59 `np_start13delay12_nc384_t` |
|---|---|---|
| start | per-unit mode and concentration (`start_dist` 3) | one Beta shared (`start_dist` 1) |
| window | [-6, 13] h (`start_below` 6) | [0, 13] h |
| ages below 0 | the cycle's end: late schizonts bursting within 6 h | none |
| delay | none | per unit, flat on [0, 12] h (`start_delay`, `delay_max` 12) |
| ring-only premise | broken on purpose | kept |
| new parameters vs 57 / 55 | none | 14 delays |

Everything else is 55-57's: no_pool, `n_c` 384, `sd_log10_total0` 1, cap 3,
whole inoculum circulating at t = 0, adapt_delta 0.95, max_treedepth 12,
seed 538065874. 59 shares the start so that the delay alone carries each
unit's phase: a per-unit mode and a per-unit delay would trade off.

**How.** Package commit cda47a0 on the scratch clone's `start-dist` branch,
installed in the scratch library only. New data `start_below`,
`start_delay`, `delay_max` in all four `archer_fit*` programs, with new
functions (`windowed_beta_starts_below`, `make_DC_gen`,
`traj_combo_partial_sum_gen`, ...) reached only when either is on, so
entries 54-57 run the code they ran before.
- An age of -x h is the compartments in the cycle's last x hours. Those
  cells are circulating at t = 0 too, which is not biology for late
  schizonts but makes no difference at the observation times (7c below).
- A delay reads the trajectory at `ts - delay` (no mortality during the
  pause; `mu` is 0). The Erlang-window truncation bounds are computed from
  `ts - delay_max`, so they cover every delay the prior allows.
- `run_check` with a delay is impractically slow (`gotchas.md`); the
  equivalent comparison is section 7f.

**Validation** (`_scripts/start-dist-validate.R` sections 2-7,
`_data/start-dist-validate-below-delay-2026-11-20.txt`):
- unchanged paths: `start_dist` 0 against the pre-start_dist build,
  difference 0; 1 against 13db2bd, 1.5e-11 on log densities of ~5e4;
  2 and 3 against fc7db77, 0 and 7.3e-12; gradients within 5e-16 relative;
- below-zero starting weights against exact pbeta, worst 1.4e-5 in
  compartment weight (shapes between 1 and 2; exact at integer shapes);
  with no below-zero part identical to the old weights (1e-16); continuous
  as cycle_length moves the lower edge across a compartment (3e-8 for a
  1e-7 relative step);
- all circulating against split by y(age) at t = 0, with half the start
  below age 0: circulating abundance at 72-216 h identical to 1.4e-15;
- layouts: 58 has 99 parameters (as 57), 59 has 87 (55's 73 plus 14);
  finite log density and gradient at 8 points; gradients against central
  differences within 1.7e-7; delay to 0 recovers 55's density, linearly
  (1.1e-3 at 1e-6 h, 1.1e-6 at 1e-9 h); `log_lik` in generated quantities
  moves exactly with the model block (differences 0 and 1.5e-11);
- fast path against a matrix exponential at the off-grid read times, 5
  layouts x `n_c` 96 and 384 x cycle lengths 35.2, 42.84 and 49.8 h,
  delays 0 to 11.5 h: worst relative difference 2.5e-12;
- `archer_stan_data()` builds both entries' data as entry 54's plus the
  start fields, `wockner-fit.R`'s helpers round-trip both, and bad inputs
  are refused.

**Runs.** SLURM 35290, array 58-59, submitted 2026-11-20 with
`WOCKFIT_LIB=/home2/lan68/plasmofit/.Rlib-dev`. Reader:
`_scripts/read-entry5859.R`.
