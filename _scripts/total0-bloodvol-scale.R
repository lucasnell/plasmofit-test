## How big is the between-individual initial-density variation that the model
## has no term for?
##
## `para` is a concentration (iRBC/mL) but `inoc_size` is a COUNT, identical
## for everyone in a cohort. The conversion is a division by blood volume,
## which differs between people and is not in the data. archer_stan_data()
## uses blood_volume_ml = 5000 for everyone, and log10_total0 is a vector over
## the 14 grp_init groups, not over the 177 series -- so there is no
## per-individual initial-density parameter at all.
##
## A CONSTANT error in blood volume is absorbed by delta_total0, as the
## package documents. Between-individual VARIATION is not: a single offset
## cannot absorb a spread. This sizes that spread against the observation
## error it would have to be detected over.
##
##   cd /home2/lan68/plasmofit/plasmofit-test
##   Rscript --vanilla _scripts/total0-bloodvol-scale.R

suppressPackageStartupMessages({library(rstan); library(readr); library(dplyr)})

## Blood volume in screened healthy adults. Nadler's equation gives roughly
## 70 mL/kg in men and 65 in women, and blood volume scales sublinearly with
## weight because adipose tissue is poorly vascularised, so CV(volume) is
## somewhat below CV(weight). A volunteer-infection cohort is weight-screened,
## so take CV in the 12-20% range and report the span rather than one number.
CV_BV <- c(0.12, 0.15, 0.20)

cat("=== between-individual spread in log10 initial density ===\n")
cat("Cells: sd on the log10 scale implied by a lognormal blood volume with\n",
    "the given coefficient of variation; sd_log10 = CV / ln(10).\n\n", sep = "")
bv <- tibble(CV_blood_volume = CV_BV, sd_log10_total0 = CV_BV / log(10))
print(as.data.frame(bv), digits = 3)

f <- readRDS("_data/wock-fit-np_wide_both.rds")
s <- summary(f)$summary
sd_i <- s[grep("^sd_iRBC\\[", rownames(s)), "mean"]

d <- read_csv("_data/wockner-cleaned.csv", col_types = "cccdcdd")
obs_per_series <- d |> count(id) |> pull(n)

cat("\n=== what it would have to be seen over ===\n")
cat(sprintf("observation error sd_iRBC: %d groups, mean %.3f, range %.3f-%.3f (log10)\n",
            length(sd_i), mean(sd_i), min(sd_i), max(sd_i)))
cat(sprintf("observations per series: median %d, range %d-%d\n",
            median(obs_per_series), min(obs_per_series), max(obs_per_series)))
se_level <- mean(sd_i) / sqrt(median(obs_per_series))
cat(sprintf("=> standard error on ONE series' mean level: %.3f log10\n", se_level))

cat("\n=== comparison ===\n")
for (cv in CV_BV) {
    sdl <- cv / log(10)
    cat(sprintf("CV %.0f%%: sd %.3f log10 = %.2f x the per-series SE, %.1f%% of observation variance\n",
                100 * cv, sdl, sdl / se_level, 100 * sdl^2 / mean(sd_i)^2))
}

cat("\nA per-series random effect on log10_total0 is the right structural fix,\n",
    "but it has to be identified against the observation error above. If the\n",
    "implied sd is well under the per-series SE, the effect is real and\n",
    "simply not estimable from these data, and adding the term would move\n",
    "variance between sd_iRBC and the new parameter without changing fit.\n", sep = "")
