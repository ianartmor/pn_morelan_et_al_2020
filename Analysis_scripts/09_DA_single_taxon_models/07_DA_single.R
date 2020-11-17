library(corncob)
library(tidyverse)
library(patchwork)
source("Analysis_scripts/00_setup.R")
source("Analysis_scripts/7_da_functions.R")


#### 16S
merge_gen_16s <- taxa_level(must_16s_merge, "Genus")
toptax_16s <- prune_taxa(top_taxa(merge_gen_16s, 12), merge_gen_16s)

### VINTAGE 16S
vintage_res_16s <- vintage_dif(toptax_16s)
vintage_plotlist_16s <- purrr::map(.x = vintage_res_16s$significant_taxa, ~ vintage_dif_single(ps = toptax_16s, ASV = .x))
names(vintage_plotlist_16s) <- vintage_res_16s$significant_taxa
vintage_plotlist_16s %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/16S/Vintage/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )

### PRECIPITATION 16S
precip_res_16s <- must_precip_growing_dif(toptax_16s)
precip_plotlist_16s <- purrr::map(.x = precip_res_16s$significant_taxa, ~ precip_dif_single(ps = toptax_16s, ASV = .x))
names(precip_plotlist_16s) <- precip_res_16s$significant_taxa
precip_plotlist_16s %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/16S/Precipitation/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )


### TA 16S
ta_res_16s <- must_ta_dif(toptax_16s)
ta_plotlist_16s <- purrr::map(.x = ta_res_16s$significant_taxa, ~ ta_dif_single(ps = toptax_16s, ASV = .x))
names(ta_plotlist_16s) <- ta_res_16s$significant_taxa
ta_plotlist_16s %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/16S/Must_TA/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )


### pH 16S
pH_res_16s <- must_pH_dif(toptax_16s)
pH_plotlist_16s <- purrr::map(.x = pH_res_16s$significant_taxa, ~ pH_dif_single(ps = toptax_16s, ASV = .x))
names(pH_plotlist_16s) <- pH_res_16s$significant_taxa
pH_plotlist_16s %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/16S/Must_pH/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )


### TSS 16S
tss_res_16s <- must_tss_dif(toptax_16s)
tss_plotlist_16s <- purrr::map(.x = tss_res_16s$significant_taxa, ~ tss_dif_single(ps = toptax_16s, ASV = .x))
names(tss_plotlist_16s) <- tss_res_16s$significant_taxa
tss_plotlist_16s %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/16S/Must_TSS/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )


#### ITS
# extract tax table
taxtab <- must_its_merge %>%
  tax_table() %>%
  data.frame()

# replace "g__"
genus_names <- taxtab$Genus %>%
  as.character() %>%
  gsub("g__", "", .)

# update tax table
taxtab$Genus <- genus_names
tax_table(must_its_merge) <- tax_table(taxtab %>% as.matrix())
merge_gen_its <- taxa_level(must_its_merge, "Genus")

toptax_its <- prune_taxa(top_taxa(merge_gen_its, 8), merge_gen_its)

### VINTAGE its
vintage_res_its <- vintage_dif(toptax_its)
vintage_plotlist_its <- purrr::map(.x = vintage_res_its$significant_taxa, ~ vintage_dif_single(ps = toptax_its, ASV = .x))
names(vintage_plotlist_its) <- vintage_res_its$significant_taxa
vintage_plotlist_its %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/ITS/Vintage/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )

### PRECIPITATION its

precip_res_its <- must_precip_growing_dif(toptax_its)
precip_plotlist_its <- purrr::map(.x = precip_res_its$significant_taxa, ~ precip_dif_single(ps = toptax_its, ASV = .x))
names(precip_plotlist_its) <- precip_res_its$significant_taxa
precip_plotlist_its %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/ITS/Precipitation/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )

### TA its
# no significant taxa
# ta_res_its <- must_ta_dif(toptax_its)
# ta_plotlist_its <- purrr::map(.x = ta_res_its$significant_taxa, ~ ta_dif_single(ps = toptax_its, ASV = .x))
# names(ta_plotlist_its) <- ta_res_its$significant_taxa
# ta_plotlist_its %>%
#  map2(
#    .x = ., .y = names(.),
#    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/ITS/Must_TA/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
#  )

### pH its
pH_res_its <- must_pH_dif(toptax_its)
pH_plotlist_its <- purrr::map(.x = pH_res_its$significant_taxa, ~ pH_dif_single(ps = toptax_its, ASV = .x))
names(pH_plotlist_its) <- pH_res_its$significant_taxa
pH_plotlist_its %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/ITS/Must_pH/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )

### TSS its
tss_res_its <- must_tss_dif(toptax_its)
tss_plotlist_its <- purrr::map(.x = tss_res_its$significant_taxa, ~ tss_dif_single(ps = toptax_its, ASV = .x))
names(tss_plotlist_its) <- tss_res_its$significant_taxa
tss_plotlist_its %>%
  map2(
    .x = ., .y = names(.),
    ~ ggsave(plot = .x, filename = paste("DA_single_taxon_models/ITS/Must_TSS/", .y, ".pdf", sep = ""), width = 6.5, height = 2, units = "in")
  )

