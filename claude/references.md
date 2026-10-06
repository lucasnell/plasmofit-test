
## Greischar & Childs, "Extraordinary parasite multiplication rates in human
## malaria infections", Trends in Parasitology (opinion)

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
