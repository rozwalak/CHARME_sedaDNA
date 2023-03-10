install.packages("rioja")
install.packages("analogue")
install.packages("janitor")
install.packages("svglite")
library(vegan)
library(rioja)
library(readxl)
library(writexl)
library(funrar)
library(tidyverse)
library(analogue)
library(janitor)
library(svglite)

#Data preparation

#raw_data
df <- read_excel("input/ATP23-002_Analiza.xlsx", sheet = "Analiza V9 - combine")

#metadata
meta <- read_excel("input/ATP23-002_Analiza.xlsx", sheet = "Metadata")
ages <- read_excel("input/Dating.xlsx")
metadata <- merge(meta, ages) %>%
  select(DNA_ID, Age, Site)

#taxonomy dataframe
taxo <- df %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  select(1:8)

#only numbers about abundance
matrix <- df %>% select(1, 3:51) %>%
  column_to_rownames(var="#ASV ID")

#relative abundance of matrix
matrix_t <- t(matrix)
rel_abundance <- make_relative(matrix_t)
rel_abundance <- t(rel_abundance)
rel_abundance <- data.frame(rel_abundance) %>%
  mutate(Total = rowSums(.)) %>%
  rownames_to_column("#ASV ID") %>%
  merge(taxo)

rel_abundance_phylum_stat <- rel_abundance %>% 
  group_by(Phylum) %>% 
  summarise(Frequency = sum(Total))

################################################################################

Marian_rel_abundance_phylum <- rel_abundance %>%
  select(53, 2:51) %>%
  select(Phylum,Total,e013,e051,e052,e053,e054,e084,e085,e086,e087,e088) %>%
  drop_na() %>%
  group_by(Phylum) %>% 
  summarize_at(vars(Total,e013,e051,e052,e053,e054,e084,e085,e086,e087,e088), funs(sum))

Marian<- data.frame(t(Marian_rel_abundance_phylum)) %>%
  row_to_names(row_number = 1) %>%
  rownames_to_column("DNA_ID")
Marian <- Marian[-1, ] %>%
  merge(metadata)
Marian[2:49] <- sapply(Marian[2:49],as.numeric)

################################################################################

Borgen_rel_abundance_phylum <- rel_abundance %>%
  select(53, 2:51) %>%
  select(Phylum,Total,e028,e057,e058,e059,e082,e083) %>%
  drop_na() %>%
  group_by(Phylum) %>% 
  summarize_at(vars(Total,e028,e057,e058,e059,e082,e083), funs(sum))

Borgen <- data.frame(t(Borgen_rel_abundance_phylum)) %>%
  row_to_names(row_number = 1) %>%
  rownames_to_column("DNA_ID")
Borgen <- Borgen[-1, ] %>%
  merge(metadata)
Borgen[2:49] <- sapply(Borgen[2:49],as.numeric)

################################################################################

Sheldon_rel_abundance_phylum <- rel_abundance %>%
  select(53, 2:51) %>%
  select(Phylum,Total,e062,e063,e064,e065,e066,e069,e070,e071,e072,e073,e074,e075,e076,e077,e078,e079,e080,e081) %>%
  drop_na() %>%
  group_by(Phylum) %>% 
  summarize_at(vars(Total,e062,e063,e064,e065,e066,e069,e070,e071,e072,e073,e074,e075,e076,e077,e078,e079,e080,e081), funs(sum))

Sheldon <- data.frame(t(Sheldon_rel_abundance_phylum)) %>%
  row_to_names(row_number = 1) %>%
  rownames_to_column("DNA_ID")
Sheldon <- Sheldon[-1, ] %>%
  merge(metadata)
Sheldon[2:49] <- sapply(Sheldon[2:49],as.numeric)

################################################################################
#Rioja Marian

Marian.rioja <- Marian %>%
  arrange(Age) %>%
  select(where(~ any(. != 0))) %>%
  select(-DNA_ID,-Age,-Site) 

y.scale <- c(1972,
             1976,
             1984,
             1995,
             1999,
             2001,
             2005,
             2011,
             2015,
             2017)

y.tks <- c(1972:2017)

p.col <- c(rep("forestgreen", times=7), rep("gold2", times=20))
ma.dist <- vegdist(Marian.rioja, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
ma.chclust <- chclust(ma.dist, method="coniss")

pol.plot <- strat.plot(Marian.rioja, yvar=y.scale, y.tks=y.scale, y.rev=FALSE, plot.line=TRUE, plot.poly=TRUE, plot.bar=FALSE, col.bar=p.col, lwd.bar=0.001, scale.percent=FALSE, xSpace=0.005, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2, title = "Marian Cove", ylabel = "Age",
                       srt.xlabel=70, clust=ma.chclust, ylabPos=3)

Marian_rioja <- addClustZone(pol.plot, ma.chclust, nZone=2, lwd=1.5, lty=2, col="grey25")

svg("output/Marian.svg")

################################################################################

#Rioja Borgen

Borgen.rioja <- Borgen %>%
  arrange(Age) %>%
  select(where(~ any(. != 0))) %>%
  select(-DNA_ID,-Age,-Site) 

y.scale <- c(1966,
             1968,
             1969,
             1970,
             1998,
             2014)

y.tks <- c(1966:2014)

p.col <- c(rep("forestgreen", times=7), rep("gold2", times=20))
ma.dist <- vegdist(Borgen.rioja, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
ma.chclust <- chclust(ma.dist, method="coniss")

pol.plot <- strat.plot(Borgen.rioja, yvar=y.scale, y.tks=y.scale, y.rev=FALSE, plot.line=TRUE, plot.poly=TRUE, plot.bar=FALSE, col.bar=p.col, lwd.bar=0.001, scale.percent=FALSE, xSpace=0.005, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2, title = "B??rgen Bay", ylabel = "Age",
                       srt.xlabel=70, clust=ma.chclust, ylabPos=3)

Borgen_rioja <- addClustZone(pol.plot, ma.chclust, nZone=2, lwd=1.5, lty=2, col="grey25")

svg("output/Borgen.svg")



################################################################################
#Rioja Sheldon

Sheldon.rioja <- Sheldon %>%
  arrange(Age) %>%
  select(where(~ any(. != 0))) %>%
  select(-DNA_ID,-Age,-Site) 

y.scale <- c(1948,
             1950,
             1953,
             1956,
             1958,
             1961,
             1964,
             1967,
             1970,
             1972,
             1973,
             1975,
             1979,
             1986,
             1993,
             2000,
             2007,
             2014)
             
y.tks <- c(1948:2014)

p.col <- c(rep("forestgreen", times=7), rep("gold2", times=20))
ma.dist <- vegdist(Sheldon.rioja, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
ma.chclust <- chclust(ma.dist, method="coniss")

pol.plot <- strat.plot(Sheldon.rioja, yvar=y.scale, y.tks=y.scale, y.rev=FALSE, plot.line=TRUE, plot.poly=TRUE, plot.bar=FALSE, col.bar=p.col, lwd.bar=0.001, scale.percent=FALSE, xSpace=0.005, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2, title = "Sheldon Cove", ylabel = "Age",
                       srt.xlabel=70, clust=ma.chclust, ylabPos=3)

Sheldon_rioja <- addClustZone(pol.plot, ma.chclust, nZone=3, lwd=1.5, lty=2, col="grey25")

svg("output/Sheldon.svg")


################################################################################

#Ochrophyta - Sheldon

Sheldon_ochrophyta <- rel_abundance %>%
  merge(taxo) %>%
  filter(Phylum == "Ochrophyta") %>%
  select(Genus,Total,e062,e063,e064,e065,e066,e069,e070,e071,e072,e073,e074,e075,e076,e077,e078,e079,e080,e081) %>%
  drop_na() %>%
  group_by(Genus) %>% 
  summarize_at(vars(Total,e062,e063,e064,e065,e066,e069,e070,e071,e072,e073,e074,e075,e076,e077,e078,e079,e080,e081), funs(sum))

Sheldon_ochrophyta <- data.frame(t(Sheldon_ochrophyta)) %>%
  row_to_names(row_number = 1) %>%
  rownames_to_column("DNA_ID")
Sheldon_ochrophyta <- Sheldon_ochrophyta[-1, ] %>%
  merge(metadata)
Sheldon_ochrophyta[2:6] <- sapply(Sheldon_ochrophyta[2:6],as.numeric)

Sheldon_ochrophyta.rioja <- Sheldon_ochrophyta %>%
  arrange(Age) %>%
  select(where(~ any(. != 0))) %>%
  select(-DNA_ID,-Age,-Site) 

y.scale <- c(1948,
             1950,
             1953,
             1956,
             1958,
             1961,
             1964,
             1967,
             1970,
             1972,
             1973,
             1975,
             1979,
             1986,
             1993,
             2000,
             2007,
             2014)

y.tks <- c(1948:2014)

p.col <- c(rep("forestgreen", times=7), rep("gold2", times=20))
ma.dist <- vegdist(Sheldon_ochrophyta.rioja, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
ma.chclust <- chclust(ma.dist, method="coniss")

pol.plot <- strat.plot(Sheldon_ochrophyta.rioja, yvar=y.scale, y.tks=y.scale, y.rev=FALSE, plot.line=TRUE, plot.poly=TRUE, plot.bar=TRUE, col.bar=p.col, lwd.bar=0.01, scale.percent=FALSE, xSpace=0.05, x.pc.lab=TRUE, x.pc.omit0=TRUE, las=2, title = "Sheldon Cove Dinoflagellates - Genera", ylabel = "Age",
                       srt.xlabel=70, ylabPos=3)

Sheldon_rioja <- addClustZone(pol.plot, ma.chclust, nZone=3, lwd=1.5, lty=2, col="grey25")

