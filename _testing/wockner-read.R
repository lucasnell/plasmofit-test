
# Much of this is from Greischar & Childs (2023)
# https://github.com/greischarlab/extraordinaryPMRs/blob/main/Plot_Fig1_S1_S6.R

suppressPackageStartupMessages({
    library(tidyverse)
})

# Note: readxl package required below, but I'm not loading it bc I only need it
# once for the `readxl::read_xlsx` function.
# Same thing for docxtractr and the `docxtractr::read_docx` and
# `docxtractr::docx_extract_tbl` functions.



wockner_dir <- "~/Box/__Cornell/_plasmofit/Wockner-2020/"


# To add inoculation sizes
inoc_size_df <- paste0(wockner_dir, "jiz557_suppl_supplementary_table_1.docx") |>
    docxtractr::read_docx() |>
    docxtractr::docx_extract_tbl(1) |>
    select(c(1, 3, 7)) |>
    slice(-1) |>
    set_names(c("trial", "cohort", "inoc_size")) |>
    mutate(inoc_size = as.numeric(inoc_size),
           cohort = paste0("CH", cohort), # << to match paras_df
           # To remove references from table:
           trial = str_remove_all(trial, "\\ \\[.*"),
           # To match values in `paras_df`:
           trial = replace_when(trial,
                               trial == "ACT-451840" ~ "ACT451840",
                               trial == "MMV048 PIB" ~ "MMV048_PIB",
                               trial == "EFITA/OZGAMe" ~ "EFITA_OZGAM",
                               trial == "DSMOZ-2" ~ "DSMOZ_2",
                               trial == "MMV048 Part B" ~ "MMV048_PartB",
                               trial == "" ~ NA)) |>
    mutate(cohort = str_remove_all(cohort, "\\/.*"), # fix EFITA_OZGAM cohorts
           # fix DSM265 cohorts:
           cohort = str_replace_all(cohort, "1A$", "1"),
           cohort = str_replace_all(cohort, "1B$", "2")) |>
    # Fill in inoculation values for empty cells
    fill(trial)



paras_df <- paste0(wockner_dir, "jiz557_suppl_supplementary_table_2.xlsx") |>
    readxl::read_xlsx(range = "R12C1:R1142C7", col_names = TRUE) |>
    set_names(c("time", "treat_day", "log10_para", "para", "tc", "trial_cohort", "id")) |>
    select(id, time, para) |>
    mutate(time = time * 24) |>  # days to hours
    # Fix weird id strings:
    mutate(id = str_replace_all(id, "EFITA:CH.?", "EFITA") |>
               str_replace_all("CH2B:CH", "CH")) |>
    mutate(id_split = str_split(id, ":"),
           trial = map_chr(id_split, \(x) x[[1]]),
           cohort = map_chr(id_split, \(x) x[[2]]),
           subject = map_chr(id_split, \(x) x[[3]])) |>
    select(-id_split) |>
    select(id, trial:subject, everything()) |>
    mutate(inoc_size = map2_dbl(trial, cohort, \(i, c) {
        inoc_size_df$inoc_size[inoc_size_df$trial == i & inoc_size_df$cohort == c]
    })) |>
    # inoc_size always groups by cohort, but not always by trial
    select(id, trial, cohort, inoc_size, subject, everything())


rm(inoc_size_df, wockner_dir)

write_csv(paras_df, "_testing/wockner-cleaned.csv")

