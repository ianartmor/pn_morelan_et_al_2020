# core analyses

#library(microbiome)
#library(microbiomeSeq)
#library(tidyverse)
#library(phyloseq)
#library(WriteXLS)
#source("setup.R")





core_analyses <- function(physeq){
  core <- core_members(physeq %>% transform("compositional"), detection = 0.0001, prevalence = 0.9)
  core_ps <- prune_taxa(core, physeq %>% transform("compositional"))
  quants <- otu_table(core_ps)  %>% as.data.frame() %>%dplyr::summarise(across(is.numeric,quantile)) %>% t %>% as.data.frame() %>% rownames_to_column("ASV") %>% as_tibble() %>% select(ASV,`0`="V1",`25`="V2", `50`="V3",`75`="V4",`100`="V5")
  taxtab <- tax_table(core_ps)[core] %>% as.data.frame %>% rownames_to_column("ASV") %>%as_tibble()
  coretab <- left_join(taxtab, quants, by="ASV") %>% as.data.frame()
  return(coretab)
}

variable_analyses <- function(physeq){
  ps_rel <- physeq %>% transform("compositional")
  core <- core_members(ps_rel, detection = 0.0001, prevalence = 0.9)
  variable <- variable_members(ps_rel, detection = 0.01, prevalence = 0.05)
  variable_ps <- prune_taxa(setdiff(variable,core), ps_rel)
  quants <- otu_table(variable_ps)  %>% as.data.frame() %>%dplyr::summarise(across(is.numeric,quantile)) %>% t %>% as.data.frame() %>% rownames_to_column("ASV") %>% as_tibble() %>% select(ASV,`0`="V1",`25`="V2", `50`="V3",`75`="V4",`100`="V5")
  taxtab <- tax_table(variable_ps) %>% as.data.frame %>% rownames_to_column("ASV") %>%as_tibble()
  variabletab <- left_join(taxtab, quants, by="ASV") %>% as.data.frame()
  return(variabletab)
}

core_all<- core_analyses(must_16s_merge)
core_2016 <- core_analyses(must_16s_2016)
core_2017 <- core_analyses(must_16s_2017)
variable_all <- variable_analyses(must_16s_merge)
metadata <- enframe(c(
              "A: Bacterial taxonomy and abundance quantiles of ASVs with mean abundances greater than 0.01% in at least 90% of samples, calculated over the entire dataset",
              "B:  Bacterial taxonomy and abundance quantiles of ASVs with mean abundances greater than 0.01% in at least 90% of samples, calculated for 2016 samples",
              "C: Bacterial taxonomy and abundance quantiles of ASVs with mean abundances greater than 0.01% in at least 90% of samples, calculated for 2017 samples",
              "D: Bacterial taxonomy and abundance quantiles of ASVs with  mean abundances greater than 1%  in between 5% and 90% of samples"), value="Tab details", name=NULL)

WriteXLS(ExcelFileName = "Figures_and_tables_check/Table_S2_core_analyses_16S.xls", x=c("core_all","core_2016","core_2017","variable_all", "metadata"))

core_all_its<- core_analyses(must_its_merge)
core_2016_its <- core_analyses(must_its_2016)
core_2017_its <- core_analyses(must_its_2017)
variable_all_its <- variable_analyses(must_its_merge)
metadata_its <- enframe(c(
                            "A: Fungal taxonomy and abundance quantiles of ASVs with mean abundances greater than 0.01% in at least 90% of samples, calculated over the entire dataset",
                            "B: Fungal taxonomy and abundance quantiles of ASVs with mean abundances greater than 0.01% in at least 90% of samples, calculated for 2016 samples",
                            "C: Fungal taxonomy and abundance quantiles of ASVs with mean abundances greater than 0.01% in at least 90% of samples, calculated for 2017 samples",
                            "D: Fungal taxonomy and abundance quantiles of ASVs with  mean abundances greater than 1%  in between 5% and 90% of samples"), 
                        value = "Tab details", name = NULL)

WriteXLS(ExcelFileName = "Figures_and_tables_check/Table_S3_core_analyses_its.xls", x=c("core_all_its","core_2016_its","core_2017_its","variable_all_its", "metadata_its"), SheetNames = c("A core_all_its","B core_2016_its","C core_2017_its","D variable_all_its", "metadata"))
