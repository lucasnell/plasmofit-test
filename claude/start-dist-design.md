# Design: a shared, one-sided starting age distribution

*Drafted 2026-11-09 with Lucas, after `findings.md`, "The inoculum's ages".
Nothing is built. The decisions at the end are his.*

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
