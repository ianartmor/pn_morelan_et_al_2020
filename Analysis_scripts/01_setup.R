# load packages and data

#Load packages
library(dunn.test)
library(dada2)
library(vegan)
library(microbiome)
library(geodist)
library(tidyverse)
library(ggordiplots)
library(corncob)
library(ggmap)
library(ggrepel)
library(shades)
library(patchwork)
library(ggsn)
library(phyloseq)
library(WriteXLS)
library(RColorBrewer)
library(knitr)
library(lme4)

# read in list of phyloseq objects
phylist <- readRDS("phyloseq_objects.rds")

# extract phyloseq objects from list
must_16s <- phylist[[1]]
must_16s_merge <- phylist[[2]]
must_its <- phylist[[2]]
must_its_merge <- phylist[[4]]



# convert merged and unmerged sample data to tibbles for easier access
samdat <- sample_data(must_16s) %>% as_tibble()
merge_samdat <- read.csv("merge_samdat.csv", row.names = 1)
merge_samdat$AVA <- factor(merge_samdat$AVA, levels=c("Santa Barbara","Monterey", "Sonoma","Mendocino","Willamette Valley"), ordered = TRUE)
merge_samdat$Shipment <- rownames(merge_samdat)
sample_data(must_16s_merge) <- sample_data(merge_samdat)
sample_data(must_its_merge) <- sample_data(merge_samdat)


# some analyses are conducted separately for each year
must_16s_2016 <- subset_samples(must_16s_merge, Year == "2016")
must_16s_2017 <- subset_samples(must_16s_merge, Year == "2017")
must_its_2016 <- subset_samples(must_its_merge, Year == "2016")
must_its_2017 <- subset_samples(must_its_merge, Year == "2017")

### plotting
# set theme
theme_set(theme_gray())
# set palette for sites
mypal <- brightness(brewer.pal(5,"Spectral"), 0.8)
# shape values
shapevals <- c(21:25)