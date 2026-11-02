
## Greischar & Childs (2023), "Extraordinary parasite multiplication rates in
## human malaria infections", Trends in Parasitology 39(8),
## doi 10.1016/j.pt.2023.05.006

`mmcm.pdf` (article) and `mmc1.pdf` (supplement) in the repo root, untracked.
Argues that established methods return implausibly large parasite
multiplication rates, and that the inflation comes from two things: some
developmental ages are easier to sample than others, and **the distribution
of developmental ages changes over the course of infection**.

**Directly relevant to thread 5.** That second mechanism is the
desynchronisation this project found `n_c` controls, and the first is a
sampling-vs-stage-structure interaction the model represents only through
`b_shape` and the Erlang chain. Read before writing up the `n_c` result:
it is prior art on the same confound, from the same department, and the
`n_c = 96` rung losing 71.8 elpd is an instance of exactly the problem it
describes. The supplement lists code and malariatherapy data at
https://github.com/greischarlab/extraordinaryPMRs (322 patients).


## Smith et al. (2020), "An intrinsic oscillator drives the blood stage cycle
## of the malaria parasite, Plasmodium falciparum", Science 368(6492), 754-759,
## doi 10.1126/science.aba4357

PDF and supplement read 2026-11-02 (`_data/new-info/`, gitignored). Fits a
phase-oscillator model to population staging data for four *P. falciparum*
strains plus HB3. **Table S6, CoV of the cycle period at each strain's best
mean period: 3D7 2.76% (39 h), SA250 10.2% (54 h), FVO 8.17% (43 h), D6
0.227% (36 h), HB3 8.27% (46 h)**; fibroblast circadian reference 8.45%. For
3D7 at +/- 1 h on the mean period, 3.65% (38 h) and 4.20% (40 h). **3D7 has
the largest model error of the five** (1.16-1.59 against 0.168-0.647), so it
is the strain their model fits worst.

**It is NOT an anchor for `n_c`, and this matters.** Their supplement states
they "do not explicitly model cycle-over-cycle variability or within-cycle
variability of an individual parasite's cycle length ... and instead assume
each parasite has a constant natural IEC frequency", their data being under
two cycles long. So their CoV is PERSISTENT between-parasite heterogeneity,
whose population spread grows LINEARLY in cycles. The Erlang chain's
CV = 1 / sqrt(`n_c`) is the per-cycle, within-lineage kind, whose spread
grows as sqrt(cycles). Reading 2.76% as 1 / sqrt(`n_c`) would give `n_c`
1313; that reading is wrong. See `_scripts/unitsim-external-anchors.R` and
`findings.md`, "What the literature can and cannot anchor".

## Park, Huang, Walzer, You, Chi & Buchler (2020), "Heterogeneous timing of
## asexual cycles in Plasmodium falciparum quantified by extended time-lapse
## microscopy", bioRxiv 10.1101/344812; data at doi 10.5061/dryad.7sqv9s4q0

PDF read 2026-11-02 (`_data/new-info/`). **The one source that separates the
two kinds of variability, by direct single-cell time-lapse of 3D7HT-GFP**,
the trials' 3D7 background. It times each founder's burst and each of its
progeny's next burst, so the burst-to-burst interval is one parasite's cycle,
and decomposes it as tau_ik = tau_i + tau_k: **sigma_i = founder variance
(between lineages, persistent) and sigma_k = sibling variance (within one
founder's progeny, per cycle)**. Table 1: synchronized rep 1 tau 54.0 h,
sigma_i 2.9, sigma_k 1.7; unsynchronized rep 1 tau 57.5, sigma_i 4.9,
sigma_k 5.0; synchronized rep 2 tau 57.3, sigma_i 2.7, sigma_k 2.7.

**sigma_k is the quantity `n_c` represents**, giving `n_c` 1009 and 450 in
the two synchronized replicates (132 unsynchronized). **sigma_i is a
mechanism the model has no term for** and is as large as sigma_k. Caveats:
in vitro; mean periods 54-57 h under the microscope against ~39-43 h
expected, which they attribute to culture conditions, having ruled out GFP
and photo-toxicity (Fig. S4, control 48.7 +/- 5.0 h); CV is scale-free but
the stress may inflate it.

## Greischar, Savill, Reece & Mideo (2024), "How to quantify developmental
## synchrony in malaria parasites", Frontiers in Malaria,
## doi 10.3389/fmala.2024.1386266

Found 2026-11-02. Fits a Leslie matrix to iRBC time series sampled several
times per cycle, with a **symmetric Beta initial age distribution and a
shape parameter** — the same device as this project's `b_shape` — and
defines synchrony as one minus the fraction of a cycle over which 99% of
iRBCs burst. **Relevant limitations they report, which match ours**:
initial synchrony and its decay are not separated (synchrony is assumed
constant); with even sampling twice per cycle, median age and bursting
duration trade off and are not identifiable; validation is on simulated
data only, because the stage-percentage method is biased and there is no
empirical benchmark. Also cites Greischar et al. (2019), Trends Parasitol.
35, 341-355, "The challenge of quantifying synchrony in malaria parasites".


## Wöckner et al. (2020), J Infect Dis 221:963-972 — the source data

PDF and supplements read 2026-11-02 (`_data/new-info/Wockner-2020/`). Useful
facts for anchoring, found while looking for one on `sd_iRBC`:
- **No replicate spread is available.** Samples were run in duplicate or
  triplicate in a single assay and **geometrically averaged on the log10
  scale** before release; Supplementary Table 2 (the released data, with
  `log10.av.ppml` and the same `Trial:CHn:Snnn` keys as
  `_data/wockner-cleaned.csv`) holds only the averages. The qPCR technical
  error that would anchor `sd_iRBC` is therefore not in the published data.
  Limit of detection 64 parasites/mL.
- Their sine-wave fit gives a parasite life cycle of **38.8 h (95% CI
  38.3-39.2)** overall, and by subject mean 39.5 h with sd 4.04 over 128
  subjects, range 30.6-55.0 (Supplementary Table 5). That sd is **between
  hosts**, not between parasites, so it anchors neither `n_c` nor `b_shape`.
- All inocula are 3D7 from the QIMR-B master cell bank; "the inoculum was not
  synchronized between cohorts", which is about phase, not spread.
