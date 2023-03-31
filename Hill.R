library(writexl)
library(readxl)
library(tidyverse)
library(funrar)
library(hillR)
library(analogue)
library(gridExtra)
library(ggpubr)
library(cowplot)

df <- read_excel("input/ATP23-012_Analiza.xlsx", sheet = "Analiza-merged")
meta <- read_excel("input/ATP23-012_Analiza.xlsx", sheet = "Metadata")
ages <- read_excel("input/Dating_biomarkers.xlsx")
metadata <- merge(meta, ages) %>%
  select(DNA_ID, Age, Site)

#select columns from specific sites
df_Marian <- df %>% select(e013,e051,e052,e053,e054,e084,e085,e086,e087,e088)
df_Borgen <- df %>% select(e028,e057,e058,e059,e082,e083,e089,e090,e091,e092,e093,e094)
df_Sheldon <- df %>% select(e062,e063,e064,e065,e066,e069,e070,e071,e072,e073,e074,e075,e076,e077,e078,e079,e080,e081)

#transposition 
df_Marian_trans <- t(df_Marian)
df_Borgen_trans <- t(df_Borgen)
df_Sheldon_trans <- t(df_Sheldon)

#calculating Hill numbers
hill_Marian <- data.frame(hill_taxa(df_Marian_trans, q = 0)) %>%
  add_column(hill_taxa(df_Marian_trans, q = 1)) %>%
  add_column(hill_taxa(df_Marian_trans, q = 2)) %>%
  rownames_to_column("DNA_ID") %>%
  merge(metadata) %>%
  arrange(Age) %>%
  rename(Hill_N0 = 2, Hill_N1 =3, Hill_N2 =4) %>%
  select(-Site) %>%
  add_column(Site = "Marian Cove")

hill_Borgen <- data.frame(hill_taxa(df_Borgen_trans, q = 0)) %>%
  add_column(hill_taxa(df_Borgen_trans, q = 1)) %>%
  add_column(hill_taxa(df_Borgen_trans, q = 2)) %>%
  rownames_to_column("DNA_ID") %>%
  merge(metadata) %>%
  arrange(Age) %>%
  rename(Hill_N0 = 2, Hill_N1 =3, Hill_N2 =4) %>%
  select(-Site) %>%
  add_column(Site = "Börgen Bay")

hill_Sheldon <- data.frame(hill_taxa(df_Sheldon_trans, q = 0)) %>%
  add_column(hill_taxa(df_Sheldon_trans, q = 1)) %>%
  add_column(hill_taxa(df_Sheldon_trans, q = 2)) %>%
  rownames_to_column("DNA_ID") %>%
  merge(metadata) %>%
  arrange(Age) %>%
  rename(Hill_N0 = 2, Hill_N1 =3, Hill_N2 =4) %>%
  select(-Site) %>%
  add_column(Site = "Sheldon Cove")

hill_Antarctica <- rbind(hill_Marian,hill_Borgen,hill_Sheldon)

write_xlsx(hill_Antarctica, "output/hill_Antarctica.xlsx")
##Create Hill-N0 and Hill-N1 in the same plot

#hill_Antarctica_gathered <- gather(hill_Antarctica, "Hill_type", "Hill_value", -DNA_ID, -Age, -Site, -Site.1)

#Hill_plot_combined <- ggplot(hill_Antarctica_gathered, aes(x=Age, y=Hill_value, colour=Hill_type)) +
#  geom_smooth(linewidth=1.3) +
#  geom_point(size=2) +
#  facet_wrap(~factor(Site, levels=c('Marian Cove', 'Börgen Bay', 'Sheldon Cove')), scales = "free_y") +
#  labs(x = "Hill-N0", y = "Age") +
#  theme_bw()
#Hill_plot_combined

#Plotting

Hill_N0_plot <- ggplot(hill_Antarctica, aes(y = Hill_N0, x = Age)) +
  geom_smooth(size=1.2, colour='black') +
  geom_point(size=1.3) +
  facet_wrap(~factor(Site, levels=c('Marian Cove', 'Börgen Bay', 'Sheldon Cove')), scales = "fixed", ncol = 1) +
  labs(x = "Age", y = "Hill-N0") +
  theme_bw() +
  ggtitle("Richness") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
Hill_N0_plot

Hill_N1_plot <- ggplot(hill_Antarctica, aes(y = Hill_N1, x = Age)) +
  geom_smooth(size=1.2, colour='black') +
  geom_point(size=1.3) +
  facet_wrap(~factor(Site, levels=c('Marian Cove', 'Börgen Bay', 'Sheldon Cove')), scales = "fixed", ncol = 1) +
  labs(x = "Age", y = "Hill-N1") +
  theme_bw() +
  ggtitle("Shannon") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
Hill_N1_plot

Hill_N2_plot <- ggplot(hill_Antarctica, aes(y = Hill_N2, x = Age)) +
  geom_smooth(size=1.2, colour='black') +
  geom_point(size=1.3) +
  facet_wrap(~factor(Site, levels=c('Marian Cove', 'Börgen Bay', 'Sheldon Cove')), scales = "free_y", ncol = 1) +
  labs(x = "Age", y = "Hill-N2") +
  theme_bw() +
  ggtitle("Simpson") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
Hill_N2_plot

#combine two plots
Hill_plot <- ggarrange(Hill_N0_plot, Hill_N1_plot, ncol=2)

Hill_plot
#save plot file
ggsave("Hill_N0_plot.png", plot = Hill_N0_plot , bg = "white", width = 7, height = 9, device = "png", path = "output")
ggsave("Hill_N1_plot.png", plot = Hill_N1_plot , bg = "white", width = 7, height = 9, device = "png", path = "output")
ggsave("Hill_N2_plot.png", plot = Hill_N2_plot , bg = "white", width = 7, height = 9, device = "png", path = "output")
ggsave("Hill_plot.png", plot = Hill_plot, bg = "white", width = 10, height = 10, device = "png", path = "output")
