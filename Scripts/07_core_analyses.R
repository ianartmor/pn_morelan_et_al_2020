# core analyses

core_analyses <- function(physeq){
  core <- core_members(physeq %>% transform("compositional"), detection = 0.0001, prevalence = 0.9)
  core_ps <- prune_taxa(core, physeq %>% transform("compositional"))
  quants <- otu_table(core_ps)  %>% as.data.frame() %>%dplyr::summarise(across(is.numeric,quantile)) %>% t %>% as.data.frame() %>% rownames_to_column("ASV") %>% as_tibble() %>% select(ASV,`0`="V1",`25`="V2", `50`="V3",`75`="V4",`100`="V5")
  taxtab <- tax_table(core_ps)[core] %>% as.data.frame %>% rownames_to_column("ASV") %>%as_tibble()
  coretab <- left_join(taxtab, quants, by="ASV") %>% as.data.frame()
  return(coretab)
}

noncore_analyses <- function(physeq){
  ps_rel <- physeq %>% transform("compositional")
  core <- core_members(ps_rel, detection = 0.0001, prevalence = 0.9)
  noncore <- variable_members(ps_rel, detection = 0.01, prevalence = 0.05)
  noncore_ps <- prune_taxa(setdiff(noncore,core), ps_rel)
  quants <- otu_table(noncore_ps)  %>% as.data.frame() %>%dplyr::summarise(across(is.numeric,quantile)) %>% t %>% as.data.frame() %>% rownames_to_column("ASV") %>% as_tibble() %>% select(ASV,`0`="V1",`25`="V2", `50`="V3",`75`="V4",`100`="V5")
  taxtab <- tax_table(noncore_ps) %>% as.data.frame %>% rownames_to_column("ASV") %>%as_tibble()
  noncoretab <- left_join(taxtab, quants, by="ASV") %>% as.data.frame()
  return(noncoretab)
}

core_all<- core_analyses(must_16s_merge)
core_2016 <- core_analyses(must_16s_2016)
core_2017 <- core_analyses(must_16s_2017)
noncore_all <- noncore_analyses(must_16s_merge)

WriteXLS(ExcelFileName = "Figures_and_tables_check/Supplemental_table_1_core_analyses_16S.xls", x=c("core_all","core_2016","core_2017","noncore_all"))

core_all_its<- core_analyses(must_its_merge)
core_2016_its <- core_analyses(must_its_2016)
core_2017_its <- core_analyses(must_its_2017)
noncore_all_its <- noncore_analyses(must_its_merge)

WriteXLS(ExcelFileName = "Figures_and_tables_check/Supplemental_table_2_core_analyses_its.xls", x=c("core_all_its","core_2016_its","core_2017_its","noncore_all_its"))
