
dif_abund <- function(physeq,n){
  
  toptax <- prune_taxa(top_taxa(physeq, n),physeq)
  
  vintage_dif <- differentialTest(formula = ~Year+AVA,
                                  phi.formula = ~ AVA,
                                  formula_null = ~ AVA,
                                  phi.formula_null = ~ AVA,
                                  data = toptax,
                                  test = "LRT", boot = FALSE,
                                  fdr_cutoff = 0.05)
  
  must_precip_growing_dif <- differentialTest(formula = ~ Precipitation.growing+AVA+Year,
                                              phi.formula = ~ AVA+Year,
                                              formula_null = ~ AVA+Year,
                                              phi.formula_null = ~ AVA+Year,
                                              data = toptax,
                                              test = "LRT", boot = FALSE,
                                              fdr_cutoff = 0.05)
  
  must_ta_dif <- differentialTest(formula = ~Must_TA+Year+AVA,
                                  phi.formula = ~ AVA+Year,
                                  formula_null = ~ AVA+Year,
                                  phi.formula_null = ~ AVA+Year,
                                  data = toptax,
                                  test = "LRT", boot = FALSE,
                                  fdr_cutoff = 0.05)
  
  must_ph_dif <- differentialTest(formula = ~Must_pH+Year+AVA,
                                  phi.formula = ~ AVA+Year,
                                  formula_null = ~ AVA+Year,
                                  phi.formula_null = ~ AVA+Year,
                                  data = toptax,
                                  test = "LRT", boot = FALSE,
                                  fdr_cutoff = 0.05)
  
  must_tss_dif <- differentialTest(formula = ~Must_TSS+Year+AVA,
                                   phi.formula = ~ AVA+Year,
                                   formula_null = ~ AVA+Year,
                                   phi.formula_null = ~ AVA+Year,
                                   data = toptax,
                                   test = "LRT", boot = FALSE,
                                   fdr_cutoff = 0.05)
  
  
  return(list( vintage_dif, must_precip_growing_dif, must_ta_dif, must_ph_dif, must_tss_dif))
}

bact_difs <- dif_abund(must_16s_merge %>% tax_glom("Genus"), 12)

sig_vintage <- must_16s_merge%>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[1]]$significant_taxa,.)

sig_vintage %>%psmelt %>% filter(!grepl("f__", Genus)) %>% group_by(Genus, Year) %>% summarize(mean(Abundance),sd(Abundance))

dif1 <- sig_vintage %>%psmelt %>% filter(!grepl("f__", Genus)) %>%ggplot(aes(x=as.factor(Year), y=Abundance, group=Year))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~Genus,nrow = 1)+scale_y_continuous(trans="log10")+xlab("")+ylab("URA")


must_precip_sig <- must_16s_merge%>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[2]]$significant_taxa,.)
bact_difs[[2]]$significant_models

dif2 <- must_precip_sig %>%psmelt %>% filter(!grepl("f__", Genus)) %>%ggplot(aes(x=Precipitation.growing, y=Abundance))+geom_point()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+xlab("Growing Season Precipitation (mm)")+ylab("URA")

must_ta_sig <- must_16s_merge%>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[3]]$significant_taxa,.)
bact_difs[[3]]$significant_models


dif3 <- must_ta_sig %>%psmelt %>% filter(!grepl("f__", Genus)) %>%ggplot(aes(x=Must_TA, y=Abundance))+geom_point()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+xlab("TA (g/L)")+ylab("URA")


must_ph_sig <- must_16s_merge%>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[4]]$significant_taxa,.)
dif4 <- must_ph_sig %>%psmelt %>% filter(!grepl("f__", Genus)) %>%ggplot(aes(x=Must_pH, y=Abundance))+geom_point()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+xlab("pH")+ylab("URA")


must_tss_sig <- must_16s_merge%>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[5]]$significant_taxa,.)
dif5 <- must_tss_sig %>%psmelt %>% filter(!grepl("f__", Genus)) %>%ggplot(aes(x=Must_TSS, y=Abundance))+geom_point()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+xlab("TSS (°Brix)")+ylab("URA")

bact_difabund <- dif1/dif2/dif3/dif4/dif5+plot_annotation(tag_levels = "A")

ggsave("Figures_and_tables_check/Figure_6_bacteria_differential_abundance.pdf", bact_difabund, width = 8.5,height = 11, units = "in" )



taxtab <- must_its_merge %>% tax_table() %>% data.frame
genus_names <- taxtab$Genus %>% as.character() %>% gsub("g__", "", .)
taxtab$Genus <- genus_names
rownames(taxtab)

tax_table(must_its_merge) <- tax_table(taxtab %>% as.matrix)

its_difs <- dif_abund(must_its_merge %>%tax_glom("Genus"),8)

sig_vintage_its <- must_its_merge %>% subset_samples(., Shipment != "CRN1_2017")%>% tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(its_difs[[1]]$significant_taxa,.)
sig_vintage_its %>%psmelt %>% filter(!grepl("f__", Genus)) %>% group_by(Genus, Year) %>% summarize(mean(Abundance),sd(Abundance))

dif1_its <- sig_vintage_its %>%psmelt  %>%ggplot(aes(x=as.factor(Year), y=Abundance, group=Year))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+ylab("URA")+xlab("")

must_precip_sig_its <- must_its_merge%>%subset_samples(., Shipment != "CRN1_2017")%>% tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(its_difs[[2]]$significant_taxa,.)
dif2_its <- must_precip_sig_its %>%psmelt %>%ggplot(aes(x=Precipitation.growing, y=Abundance))+geom_point()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+ylab("URA")+xlab("Growing Season Precipitation (mm)")

must_ph_sig_its <- must_its_merge%>%subset_samples(., Shipment != "CRN1_2017")%>%tax_glom("Genus") %>%  transform("compositional") %>%  prune_taxa(its_difs[[4]]$significant_taxa,.)
dif4_its <- must_ph_sig_its %>%psmelt  %>%ggplot(aes(x=Must_pH, y=Abundance))+geom_point()+geom_smooth(method="lm",se=FALSE)+scale_y_continuous(trans="log10")+ylab("URA")+xlab("pH")+facet_wrap(~Genus)

must_tss_sig_its <- must_its_merge %>%subset_samples(., Shipment != "CRN1_2017")%>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(its_difs[[5]]$significant_taxa,.)
dif5_its <- must_tss_sig_its %>%psmelt %>%ggplot(aes(x=Must_TSS, y=Abundance))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+xlab("TSS (°Brix)")+ylab("URA")+facet_wrap(~Genus)

fung_difabund <- dif1_its/dif2_its/dif4_its/dif5_its+plot_annotation(tag_levels = "A")

ggsave("Figures_and_tables_check/Figure_7_fungi_differential_abundance.pdf", fung_difabund, width = 8.5,height = 11, units = "in"  )



sig_vintage_its <- must_its_merge %>% tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(its_difs[[1]]$significant_taxa,.)
sig_vintage_its %>%psmelt %>% filter(!grepl("f__", Genus)) %>% group_by(Genus, Year) %>% summarize(mean(Abundance),sd(Abundance))

dif1_its <- sig_vintage_its %>%psmelt  %>%ggplot(aes(x=as.factor(Year), y=Abundance, group=Year))+geom_point()+geom_boxplot(alpha=0.5)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+ylab("URA")+xlab("")

must_precip_sig_its <- must_its_merge%>% tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(its_difs[[2]]$significant_taxa,.)
dif2_its <- must_precip_sig_its %>%psmelt %>%ggplot(aes(x=Precipitation.growing, y=Abundance))+geom_point()+geom_smooth(method="lm", se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+ylab("URA")+xlab("Growing Season Precipitation (mm)")

must_ph_sig_its <- must_its_merge%>%tax_glom("Genus") %>%  transform("compositional") %>%  prune_taxa(its_difs[[4]]$significant_taxa,.)
dif4_its <- must_ph_sig_its %>%psmelt  %>%ggplot(aes(x=Must_pH, y=Abundance))+geom_point()+geom_smooth(method="lm",se=FALSE)+scale_y_continuous(trans="log10")+ylab("URA")+xlab("pH")+facet_wrap(~Genus)

must_tss_sig_its <- must_its_merge %>%tax_glom("Genus") %>% transform("compositional") %>%  prune_taxa(its_difs[[5]]$significant_taxa,.)
dif5_its <- must_tss_sig_its %>%psmelt %>%ggplot(aes(x=Must_TSS, y=Abundance))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~Genus)+scale_y_continuous(trans="log10")+xlab("TSS (°Brix)")+ylab("URA")+facet_wrap(~Genus)

fung_difabund <- dif1_its/dif2_its/dif4_its/dif5_its+plot_annotation(tag_levels = "A")

ggsave("Figures_and_tables_check/Supplemental_figure_3_fungi_differential_abundance_unmodified.pdf", fung_difabund, width = 8.5,height = 11, units = "in"  )
