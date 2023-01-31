devtools::install_github("HOPE-UIB-BIO/R-Ratepol-package")


library(RRatepol)
library(tidyverse)

example_data <-  
  RRatepol::example_data

dplyr::glimpse(example_data)

example_data$sample_age[[1]]

sequence_02 <-
  RRatepol::fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]],
    smooth_method = "shep",
    DC = "bray",
    Working_Units = "levels",
    standardise = TRUE,
    N_individuals = 150,
    rand = 1000,
    use_parallel = TRUE)
example_data$pollen_data[[1]]

RRatepol::fc_plot_RoC_sequence(
  data_source = sequence_02)
example_data$sample_age[[1]]


sequence_03 <-
  RRatepol::fc_estimate_RoC(
    data_source_community = example_data$pollen_data[[1]],
    data_source_age = example_data$sample_age[[1]],
    smooth_method = "shep",
    DC = "bray",
    Working_Units = "MW",
    bin_size = 500,
    Number_of_shifts  = 5,
    standardise = TRUE,
    N_individuals = 150,
    rand = 1000,
    use_parallel = TRUE)

RRatepol::fc_plot_RoC_sequence(
  data_source = sequence_03)


biomarkers <- read.csv("m2-biomarkers.csv")

age <- read.csv("age_depth.csv")
age$age <- as.numeric(age$age)


gl_roc <-
  RRatepol::fc_estimate_RoC(
    data_source_community = biomarkers,
    data_source_age = age,
    smooth_method = "shep",
    DC = "euc",
    Working_Units = "MW",
    bin_size = 5,
    Number_of_shifts  = 5,
    standardise = TRUE,
    N_individuals = 18,
    rand = randomisations,
    use_parallel = TRUE # use_parallel = TRUE to use parallel calculation
  ) 
gl_roc


install.packages("hillR")
library(hillR)

dummy = FD::dummy
comm = dummy$abun

hill <- hill_taxa(comm, q = 0)
df <- data.frame(hill, stringsAsFactors = TRUE)

hill_taxa_parti(comm, q = 0)
