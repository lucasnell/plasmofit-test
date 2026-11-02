
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

Found 2026-11-02 while answering whether any model output could choose `n_c`.
Fits a phase-oscillator model to population staging data for four *P.
falciparum* strains (plus HB3 from earlier work) and reports a **coefficient
of variation of the cycle period of 0.23%-10.18% of the mean, depending on
strain**; the circadian reference they compare with is 8.45%. Per-strain
means and CVs are in their Table S6, NOT yet read, so **3D7's value is
unknown to us** — that is the number this project needs.

**Why it matters.** In the Erlang chain a lineage crosses `n_c` exponential
compartments per cycle, so the model implies CV = 1 / sqrt(`n_c`)
(`_scripts/unitsim-external-anchors.R`). That makes `n_c` externally
checkable: `n_c` 96 implies 10.2%, 384 implies 5.1%, 1024 implies 3.1%.
Caveat: their CV is itself a model estimate from population data, not a
single-cell measurement, and in vitro.

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
