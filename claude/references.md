
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
the two synchronized replicates and **132 in the unsynchronized one**.
**sigma_i is a mechanism the model has no term for** and is as large as
sigma_k.

**Why the synchronized replicates are the right analogue** (Lucas,
2026-11-09, and it is Park's own explanation). Their Discussion: "artificial
synchronization reduced both sibling and founder variation in the asexual
cycle period by two-fold (Table 1). The mechanism for this variance
reduction is unknown, but **selection for ring-stage parasites** could
reduce asexual cycle variability by affecting cell physiology and/or the
population demographics". The Wockner inoculum is ring-selected twice over
-- venous blood holds only circulating iRBCs, and only rings survive
Glycerolyte freezing -- so it went through the step Park names.

**This ground is chosen because it is NOT circular**, and that is the whole
point of an external anchor. The tempting alternative, "IBSM infections
have a high sine-wave amplitude so they are synchronous so use the
synchronized replicates", uses the data's own synchrony to pick the anchor
that is meant to constrain synchrony. Ring selection is a fact about how the
inoculum was prepared, known independently of any fit.

Caveats, in order of how much they could matter:
1. **The two-fold difference is the finding; ring selection is Park's
   speculation about it**, offered in discussion and not tested.
2. **It may not persist, and NO AVAILABLE DATA CAN SETTLE THAT** (checked
   2026-11-09). Park's imaging "was limited to 80+ hours because uninfected
   red blood cells started lysing", against a ~55 h cycle, so each lineage
   is followed for ONE burst-to-burst interval: founder's first burst to its
   progeny's second burst, and no third burst. Their Dataset S1 (Dryad
   doi 10.5061/dryad.7sqv9s4q0) holds exactly those three columns, so
   downloading it would not help; "quantification of variability in the
   timing of bursts over one asexual cycle" is their own description. A
   within-dataset decay test is not available either: in the synchronized
   culture all founders burst within ~10 h of release (mean 6.8, sd 4.1),
   so time since ring selection barely varies. Smith et al.'s data are
   under two cycles likewise. These observations run 72-216 h, cycles 1.6
   to 4.8 after injection, so if the reduction is transient the rate would
   rise over the infection and `n_c`, one constant for the whole
   trajectory, would average over it, weighted toward the later cycles
   where the data are. The decay-law test (2026-10-08, SLURM 29684) is the
   closest the real data come, and it could not distinguish sqrt(cycles)
   from a linear component, so they lack the power to detect it either.
   **Treat this as a standing limitation, not an open task.**
3. In vitro; mean periods 54-57 h under the microscope against ~39-43 h
   expected, which they attribute to culture conditions, having ruled out
   GFP burden and photo-toxicity (Fig. S4, control 48.7 +/- 5.0 h). CV is
   scale-free, but culture stress could inflate it.
4. Three replicates, no interval on sigma_k.

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


## The QIMR-B 3D7 inoculum: McCarthy et al. (2011), PLoS One 6:e21914;
## Cheng et al. (1997), Am J Trop Med Hyg 57:495-500

Traced 2026-11-09 from Wöckner et al. (refs 9 and 37; their Discussion says
all IBSM inocula "were prepared from a single donor unit collected from a
volunteer experimentally infected by mosquito bite from a single P falciparum
3D7 in vitro culture"). McCarthy et al. 2011 (Methods, read via PMC
PMC3159571): mosquitoes fed on 3D7 gametocytes bit a healthy volunteer; **a
500 mL unit of blood was taken 13 days after the bites, six hours after he
became ill with high fever**; leucocyte-depleted, mixed with **Glycerolyte
57**, frozen in 1 mL vials in liquid nitrogen. Inocula: an aliquot thawed,
washed, diluted; **thaw to injection <= 60 min**, on ice; "~30% of blood
stage parasites were viable". Donor parasitaemia and stage composition are
NOT reported there; they are in Cheng et al. 1997, below.

**Cheng et al. 1997, read 2026-11-09** (PDF from Lucas, `_data/malaria-inoc-x/`,
gitignored; scanned, read as page images). Two donors, bitten on day 0 and
again on day 2 by sporozoite-positive mosquitoes; parasites first
PCR-detectable on day 5.5 (volunteer 2) and 6.5 (volunteer 1). 500 mL was
taken from each "6 hr after they became ill and developed high fever (13 and
14 days after the initial mosquito bites, respectively)"; the Results date
the two collections to days 12.5 and 13.5. **Parasitaemia at collection, by
microscopy: 18 per uL (volunteer 1) and 290 per uL (volunteer 2)**; the
stored blood from volunteer 2 held 400 infected cells per uL of packed
cells. Leucocyte-filtered, Glycerolyte 57, 1 mL aliquots; thawed into 12%
then 1.6% saline, washed, injected. **Stage composition is not reported.**
Their own growth model assumed "the inoculum had two broods of parasites
separated by 24 hr and that parasites older than 24 hr would be
sequestered", under a 48 h cycle, and estimated **39% of the parasites
dividing in the first 24 hr**.

## Woolley et al. (2021), "Development and evaluation of a new Plasmodium
## falciparum 3D7 blood stage malaria cell bank for use in malaria volunteer
## infection studies", Malar J 20:93, doi 10.1186/s12936-021-03627-z

Read 2026-11-09 (PDF from Lucas). Confirms that the bank behind every study
in `_data/wockner-cleaned.csv`, **3D7-V2, was produced in 1995 from Cheng et
al.'s two volunteers**, 500 mL each "six hours following the onset of fever",
and that 3D7-V2 was chosen over 3D7-V1 for its higher parasitaemia (0.1%
against under 0.01% of erythrocytes). Two numbers bear on the inoculum's
ages:
- **Inocula are dosed on RING-STAGE parasites.** The dilution targets "the
  number of viable ring-stage parasites in each 2 ml inoculum, taking into
  account ... the percentage of parasitized erythrocytes, the percentage of
  ring-stage erythrocytes, and parasite viability", and the percentage of
  ring-stage parasites is measured by film microscopy for each bank.
- **Their new bioreactor bank was 96% ring-stage.** That is a cultured,
  deliberately staged bank, not the 1995 donor blood, so it is an analogue
  rather than a measurement of 3D7-V2; the ring percentage of 3D7-V2 itself
  is not given in either paper.

Two consequences for the model (`findings.md`, "The inoculum's ages"): the
injected iRBCs are the donor's CIRCULATING iRBCs, and only rings survive
Glycerolyte freezing (ATCC product sheet 50110, read 2026-11-09: "Only young
cells (rings) can be frozen in glycerolyte medium ... because their
membranes are more robust"); and every unit received the SAME donor
material, so all 14 units share one t = 0 age distribution.

## Cao et al. (2019), "Modeling the dynamics of Plasmodium falciparum
## gametocytes in humans during malaria infection", eLife 8:e49058

Read 2026-11-09, PDF confirming the web read. Within-host model fitted to 17
volunteers given ~2800 3D7 iRBCs (Collins et al. 2018). The inoculum's ages
are **Gaussian with mean mu and sd sigma, estimated per volunteer with flat
priors U(0, 35) h and U(0, 20) h**, no biological justification given;
cycle fixed at 42 h, sequestration age fixed at 25 h. Posteriors only in
figure supplements. Prior art for the same free initial age distribution
this project uses, equally unconstrained.
