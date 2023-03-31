install.packages("remotes")
remotes::install_github("HOPE-UIB-BIO/R-Ratepol-package",force = TRUE)
install.packages("vctrs")

library(RRatepol)
library(tidyverse)
library(readxl)
library(funrar)
trace(fc_plot_RoC_sequence, edit = T)
#Data preparation

#raw_data
df <- read_excel("input/ATP23-002_Analiza.xlsx", sheet = "Analiza V9 - combine")

#metadata
meta <- read_excel("input/ATP23-002_Analiza.xlsx", sheet = "Metadata")
ages <- read_excel("input/Dating.xlsx")
metadata <- merge(meta, ages) %>%
  select(DNA_ID, Age, Site)


matrix <- df %>% select(1, 3:51) %>%
  column_to_rownames(var="#ASV ID")

#relative abundance of matrix
matrix_t <- t(matrix)
rel_abundance <- make_relative(matrix_t)

rel_abundance <- data.frame(rel_abundance) %>%
  mutate_at(vars(), ~ . * 100) %>%
  rownames_to_column("DNA_ID")
rel_abundance <- rel_abundance %>%
  inner_join(metadata) %>%
  select(where(~ any(. != 0)))

################################################################################

Marian_rel_abundance <- rel_abundance %>%
  filter(Site == "Marian Cove") %>%
  select(where(~ any(. != 0))) %>%
  select(-Age, -Site) %>%
  rename(sample_id = DNA_ID)

Marian_age <- rel_abundance %>%
  filter(Site == "Marian Cove") %>%
  select(DNA_ID, Age) %>%
  rename(sample_id = DNA_ID) %>%
  rename(age = Age)

################################################################################
set.seed(123)
Marian_ROC <-
  RRatepol::fc_estimate_RoC(
    data_source_community = Marian_rel_abundance,
    data_source_age = Marian_age,
    smooth_method = "none",
    DC = "bray",
    Working_Units = "MW",
    bin_size = 4,
    Number_of_shifts  = 2,
    standardise = FALSE,
    rand = 1000,
    use_parallel = TRUE)
Marian_ROC

Marian_with_peaks <-
  RRatepol::fc_detect_peak_points(
    data_source = Marian_ROC,
    sel_method = "trend_non_linear")
Marian_with_peaks

fc_plot_RoC_sequence(
  Marian_ROC,
  age_threshold = NULL,
  Roc_threshold = NULL,
  Peaks = FALSE,
  trend = NULL
)

ggplot(Marian_with_peaks, aes(y = ROC, x = Age)) +
  geom_smooth(linewidth=1.6, colour='black') +
  labs(x = "Age", y = "Rate of Change") +
  theme_bw()

################################################################################

Borgen_rel_abundance <- rel_abundance %>%
  filter(Site == "Börgen Bay") %>%
  select(-Age, -Site) %>%
  rename(sample_id = DNA_ID) %>%
  select(where(~ any(. != 0)))

Borgen_age <- rel_abundance %>%
  filter(Site == "Börgen Bay") %>%
  select(DNA_ID, Age) %>%
  rename(sample_id = DNA_ID) %>%
  rename(age = Age)

################################################################################

Borgen_ROC <-
  RRatepol::fc_estimate_RoC(
    data_source_community = Borgen_rel_abundance,
    data_source_age = Borgen_age,
    smooth_method = "shep",
    DC = "bray",
    Working_Units = "MW",
    bin_size = 5,
    Number_of_shifts  = 2,
    standardise = FALSE,
    rand = 1000,
    use_parallel = TRUE)
Borgen_ROC

Borgen_with_peaks <-
  RRatepol::fc_detect_peak_points(
    data_source = Borgen_ROC,
    sel_method = "")
Borgen_with_peaks

ggplot(Borgen_with_peaks, aes(y = ROC, x = Age)) +
  geom_smooth(size=1.6, colour='black') +
  labs(x = "Age", y = "Rate of Change") +
  theme_bw()

################################################################################

Sheldon_rel_abundance <- rel_abundance %>%
  filter(Site == "Sheldon Cove") %>%
  select(-Age, -Site) %>%
  rename(sample_id = DNA_ID) %>%
  select(where(~ any(. != 0)))

Sheldon_age <- rel_abundance %>%
  filter(Site == "Sheldon Cove") %>%
  select(DNA_ID, Age) %>%
  rename(sample_id = DNA_ID) %>%
  rename(age = Age)

################################################################################
set.seed(123)
Sheldon_ROC <-
  RRatepol::fc_estimate_RoC(
    data_source_community = Sheldon_rel_abundance,
    data_source_age = Sheldon_age,
    smooth_method = "none",
    DC = "bray",
    Working_Units = "MW",
    bin_size = 3,
    Number_of_shifts  = 1,
    standardise = FALSE,
    rand = 1000,
    use_parallel = TRUE)
Sheldon_ROC

Sheldon_with_peaks <-
  RRatepol::fc_detect_peak_points(
    data_source = Sheldon_ROC,
    sel_method = "trend_non_linear")
Sheldon_with_peaks

fc_plot_RoC_sequence(
  Sheldon_with_peaks,
  age_threshold = NULL,
  Roc_threshold = NULL,
  Peaks = TRUE,
  trend = "trend_non_linear"
)

ggplot(Sheldon_with_peaks, aes(y = ROC, x = Age)) +
  geom_smooth(size=1.6, colour='black') +
  #geom_point(size=2.75) +
  labs(x = "Age", y = "Rate of Change") +
  theme_bw()

