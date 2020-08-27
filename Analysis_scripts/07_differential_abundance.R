#library(dada2)
#library(microbiomeSeq)
#library(microbiome)
#library(corncob)
#library(patchwork)
#library(cowplot)
#library(tidyverse)
#library(phyloseq)
#source("00_setup.R")



# differential abundance using corncob
dif_abund <- function(physeq,n){
  set.seed(5)
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


bact_dif <- dif_abund(must_16s_merge %>% taxa_level("Genus"), 12)
#subset to include only taxa associated with vintage
sig_vintage <- must_16s_merge%>% taxa_level("Genus") %>%transform("compositional")%>%  prune_taxa(bact_difs[[1]]$significant_taxa,.)

#sig_vintage %>%psmelt %>% filter(!grepl("f__", OTU)) %>% group_by(OTU, Year) %>% summarize(mean(Abundance),sd(Abundance))

# split data by OTU, then plot
dif1 <- sig_vintage %>% psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  purrr::map2(.x = ., .y=names(.), .f = ~ ggplot(.,aes(x=as.factor(Year), y=Abundance*100, group=Year, fill=AVA,shape=as.factor(Year)))+
                geom_point(size=2)+geom_boxplot(alpha=0.2)+
                xlab("")+ylab("")+ggtitle(.y)+
                theme(plot.title  = element_text(size = 7, hjust = -0.1),axis.title = element_text(size=9,face="plain"), legend.position = "none", 
                      panel.border = element_blank(),title = element_text(face="italic"))+ 
                scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals)+
                scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
              
  )

# add URA label to only leftmost plot
dif1[[1]] <- dif1[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))


# use modified gg.gap() function (see setup.R for code) to make gapped plots so outliers can be shown on same graph.
# rel heights parameter used to try to get ticks to be approximately the same spacing. 
dif1[[1]] <- gg.gap( dif1[[1]], ylim = c(0,16), segments = c(4,12), tick_width = 2, rel_heights = c(1,0,0.9))
dif1[[1]]<- ggdraw(dif1[[1]])+draw_label("//", x=0.52,y=0.54, angle=305)

dif1[[2]] <- gg.gap( dif1[[2]], ylim = c(0,19), segments = c(5,18), rel_heights = c(2.8,0,1), tick_width = 1,margin = c(top = 0, right = 0, bottom = 0, left = 0))
dif1[[2]]<- ggdraw(dif1[[2]])+draw_label("//", x=0.33,y=0.74, angle=305)

# rel width needs to be adjusted slightly due to leftmost plot having a y axis label
dif1 <- plot_grid(plotlist = dif1, nrow=1, rel_widths = c(1.2,1,1,1,1,1))


#extract beta binomial regression coefficients
vintage_coefficients <- bact_difs[[1]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Year")) %>% 
  mutate(rowname=bact_difs[[1]]$significant_taxa) 

#subset to include only taxa associated with precip
must_precip_sig <- must_16s_merge%>%taxa_level("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[2]]$significant_taxa,.)

# split data by OTU, then plot
dif2 <- must_precip_sig %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  purrr::map2 (.x=., .y= names(.),.f= ~ ggplot(.,aes(x=Precipitation.growing, y=Abundance*100,fill=AVA, shape=as.factor(Year)),show.legend = FALSE)+
                 geom_point()+
                 #geom_smooth(method="lm", se=FALSE)+
                 ggtitle({.y})+
                 xlab("Precipitation (mm)")+ylab("")+
                 theme(plot.title  = element_text(size = 9),axis.title = element_text(size=9,face="plain"), 
                       panel.border = element_blank(),
                       legend.position = "none",title = element_text(face="italic"))+ 
                 scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals)+
                 scale_y_continuous(labels = scales::number_format(accuracy = 0.01)))

# add URA label
dif2[[1]] <- dif2[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))

dif2 <- plot_grid(plotlist = dif2, nrow=1)

#extract beta binomial regression coefficients
precipitation_coefficients <- bact_difs[[2]]$significant_models %>%
  purrr::map(coef) %>% 
  map(data.frame) %>%
  map(rownames_to_column) %>% 
  purrr::map_df(.,.f=~filter(.,rowname=="mu.Precipitation.growing")) %>% 
  mutate(rowname=bact_difs[[2]]$significant_taxa)

# subset
must_ta_sig <- must_16s_merge%>%taxa_level("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[3]]$significant_taxa,.)

# extract coefficients
ta_coefficients <- bact_difs[[3]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Must_TA")) %>% 
  mutate(rowname=bact_difs[[3]]$significant_taxa)

# split and plot
dif3 <- must_ta_sig %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  map2(.x=.,.y=names(.),.f=~ ggplot(.,aes(x=Must_TA, y=Abundance*100,fill=AVA, shape=as.factor(Year)))+
         geom_point()+
         #geom_smooth(method="lm", se=FALSE)+
         ggtitle({.y})+
         xlab(expression(TA~(g~tartaric~acid~L^{"-1"})))+
         ylab("")+
         theme(plot.title  = element_text(size = 9),axis.title = element_text(size=9,face="plain"), 
               legend.title = element_text(size=9, face = "plain"),
               legend.text = element_text(size=7), 
               panel.border = element_blank(),
               legend.spacing.y = unit(1, "mm"), legend.box = "horizontal",title = element_text(face="italic") )+ 
         scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals, name="Vintage")+
         scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
         guides(fill = guide_legend(override.aes=list(shape=21)))
  )

# add y label
dif3[[1]] <- dif3[[1]]+ ylab("URA (%)")+theme(axis.title = element_text(size=9))

dif3 <- plot_grid(plotlist=dif3, nrow=1)

#subset
must_ph_sig <- must_16s_merge%>%taxa_level("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[4]]$significant_taxa,.)

#extract coefficients
ph_coefficients <- bact_difs[[4]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Must_pH")) %>% 
  mutate(rowname=bact_difs[[4]]$significant_taxa)

# split and plot
dif4 <- must_ph_sig %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  map2(.x=.,.y=names(.),.f=~ ggplot(.,aes(x=Must_pH, y=Abundance*100, fill=AVA, shape=as.factor(Year)))+
         geom_point()+
         #geom_smooth(method="lm", se=FALSE)+
         xlab("pH")+ylab("")+ggtitle({.y})+
         theme(plot.title  = element_text(size = 9),axis.title = element_text(size=9, face = "plain"), 
               panel.border = element_blank(),
               legend.position = "none",title = element_text(face="italic"))+ 
         scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals) +
         scale_y_continuous(labels = scales::number_format(accuracy = 0.01)))
# add y axis label
dif4[[1]] <- dif4[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))

# add gap for outlier
dif4[[1]]<- gg.gap(dif4[[1]], ylim=c(0,22.25), segments = c(2,21.75), rel_heights = c(2.4,0,1), tick_width = 0.5)
dif4[[1]]<- ggdraw(dif4[[1]])+draw_label("//", x=0.3,y=0.71, angle=305)


dif4 <- plot_grid(plotlist=dif4, nrow=1)

#subset
must_tss_sig <- must_16s_merge%>%taxa_level("Genus") %>% transform("compositional") %>%  prune_taxa(bact_difs[[5]]$significant_taxa,.)

#extract coefficient
tss_coefficients <- bact_difs[[5]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Must_TSS")) %>% 
  mutate(rowname=bact_difs[[5]]$significant_taxa)

# split and plot
dif5 <- must_tss_sig %>%psmelt %>% filter(!grepl("f__", OTU)) %>%split(.$OTU) %>% 
  map2(.x=.,.y=names(.),.f=~ggplot(.,aes(x=Must_TSS, y=Abundance*100, fill=AVA, shape=as.factor(Year)))+geom_point()+
         #geom_smooth(method="lm", se=FALSE)+
         xlab("TSS (°Brix)")+ylab("")+ggtitle({.y})+
         theme(plot.title  = element_text(size = 9),axis.title = element_text(size=9,face="plain"), 
               panel.border = element_blank(),
               legend.position = "none",title = element_text(face="italic"))+ 
         scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
         scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals))
# add y label
dif5[[1]] <- dif5[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))

# add gap for outlier
dif5[[1]]<- gg.gap(dif5[[1]], ylim=c(0,22.25), segments = c(2,21.75), rel_heights = c(2.4,0,1), tick_width = 0.5)
dif5[[1]]<- ggdraw(dif5[[1]])+draw_label("//", x=0.2,y=0.71, angle=305)

dif5 <- plot_grid(plotlist = dif5, nrow=1)

# combine plots using patchwork
bact_difabund <- dif1/dif2/dif3/dif4/dif5+plot_annotation(tag_levels = "A")+theme(axis.title = element_text(size=9))

ggsave("Figures_and_tables_check/Figure_6_bacteria_differential_abundance.pdf", bact_difabund, width = 8.5,height = 11, units = "in" )

#WriteXLS(x=c("vintage_coefficients", "precipitation_coefficients", "ta_coefficients", "ph_coefficients","tss_coefficients"), ExcelFileName = "Table_S5.xls")


# extract tax table
taxtab <- must_its_merge %>% tax_table() %>% data.frame

#replace "g__"
genus_names <- taxtab$Genus %>% as.character() %>% gsub("g__", "", .)

#update tax table
taxtab$Genus <- genus_names
tax_table(must_its_merge) <- tax_table(taxtab %>% as.matrix)

#differntial abundance
its_difs <- dif_abund(must_its_merge %>% taxa_level("Genus"), 8)

#subset
sig_vintage_its <- must_its_merge %>% taxa_level("Genus") %>%
  transform("compositional") %>%  prune_taxa(its_difs[[1]]$significant_taxa,.)

#sig_vintage_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>% group_by(OTU, Year) %>% summarize(mean(Abundance),sd(Abundance))

#extract coefficients
vintage_coefficients_its <- its_difs[[1]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Year")) %>% 
  mutate(rowname=its_difs[[1]]$significant_taxa)

#subset and plot
dif1_its <- sig_vintage_its %>% psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  purrr::map2(.x = ., .y=names(.), .f = ~ ggplot(.,aes(x=as.factor(Year), y=Abundance*100, group=Year, fill=AVA, shape=as.factor(Year)))+
                geom_point(show.legend = TRUE)+geom_boxplot(alpha=0.2, show.legend = FALSE)+
                xlab("")+ylab("")+ggtitle({.y})+
                theme(plot.title  = element_text(size = 9, face="italic"),
                      axis.title = element_text(size=9, face="plain"), 
                      panel.border = element_blank(),
                      legend.title = element_text(size=9),
                      legend.text = element_text(size=7), legend.spacing.y = unit(1, "mm"), legend.box = "horizontal" )+ 
                scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals, name="Vintage")+
                guides(fill = guide_legend(override.aes=list(shape=21)))+ 
                scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
                scale_fill_manual(values = mypal))
# add y label
dif1_its[[1]] <- dif1_its[[1]] + ylab("URA (%)")
dif1_its <- plot_grid(plotlist = dif1_its, nrow=1)

# subset
must_precip_sig_its <- must_its_merge%>% taxa_level("Genus") %>% 
  transform("compositional") %>%  prune_taxa(its_difs[[2]]$significant_taxa,.)

#split and plot
dif2_its <- must_precip_sig_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  purrr::map2 (.x=., .y= names(.),.f= ~ ggplot(.,aes(x=Precipitation.growing, y=Abundance*100, fill=AVA, shape=as.factor(Year)))+
                 geom_point()+
                 #geom_smooth(method="lm", se=FALSE)+
                 ggtitle({.y})+
                 xlab("Precipitation (mm)")+ylab("")+
                 theme(plot.title  = element_text(size = 9, face="italic"),
                       panel.border = element_blank(),
                       axis.title = element_text(size=9, face="plain"), legend.position = "none")+ 
                 scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
                 scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals))
# add y label
dif2_its[[1]] <-  dif2_its[[1]]+ylab("URA (%)") 

# add gaps for outliers
dif2_its[[2]] <-  gg.gap( dif2_its[[2]], ylim = c(0,94), segments = c(10,92), tick_width = c(2,2), rel_heights = c(3,0,1))
dif2_its[[2]]<- ggdraw(dif2_its[[2]])+draw_label("//", x=0.31,y=0.755, angle=305)

dif2_its[[3]] <- gg.gap( dif2_its[[3]], ylim = c(4,100), segments = c(6,90), tick_width = c(2,2),rel_heights = c(1,0,2.5))
dif2_its[[3]]<- ggdraw(dif2_its[[3]])+draw_label("//", x=0.335,y=0.295, angle=305)

dif2_its <- plot_grid(plotlist = dif2_its, nrow=1)

#subset
must_ph_sig_its <- must_its_merge%>%taxa_level("Genus") %>% 
  transform("compositional") %>%  prune_taxa(its_difs[[4]]$significant_taxa,.)

#extract coefficients
precipitation_coefficients_its <- its_difs[[2]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,  .f= ~filter(.,rowname=="mu.Precipitation.growing")) %>% 
  mutate(rowname=its_difs[[2]]$significant_taxa)

# split and plot
dif4_its <- must_ph_sig_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  map2(.x=.,.y=names(.),.f=~ ggplot(.,aes(x=Must_pH, y=Abundance*100, fill=AVA, shape=as.factor(Year)))+
         geom_point()+
         #geom_smooth(method="lm", se=FALSE)+
         xlab("pH")+ylab("")+ggtitle({.y})+
         theme(plot.title  = element_text(size = 9, face="italic"),axis.title = element_text(size=9, face="plain"),
               panel.border = element_blank(),
               legend.position = "none")+ 
         scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
         scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals))

# add y label
dif4_its[[1]] <-  dif4_its[[1]]+ylab("URA (%)") 

#add gap for outlier
dif4_its[[1]] <- gg.gap( dif4_its[[1]], ylim = c(4,100), segments = c(6,90), tick_width = c(2,2),rel_heights = c(1,0,2.5))
dif4_its[[1]]<- ggdraw(dif4_its[[1]])+draw_label("//", x=0.43,y=0.295, angle=305)

dif4_its <- plot_grid(plotlist=dif4_its, nrow=1)

# extract coefficients
ph_coefficients_its <- its_difs[[4]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Must_pH")) %>% 
  mutate(rowname=its_difs[[4]]$significant_taxa)

#subset
must_tss_sig_its <- must_its_merge %>%taxa_level("Genus") %>% 
  transform("compositional") %>%  prune_taxa(its_difs[[5]]$significant_taxa,.)

#split and plot
dif5_its <- must_tss_sig_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>%split(.$OTU) %>% 
  map2(.x=.,.y=names(.),.f=~ggplot(.,aes(x=Must_TSS, y=Abundance*100, fill=AVA, shape=as.factor(Year)))+geom_point()+
         # geom_smooth(method="lm", se=FALSE)+
         xlab("TSS (°Brix)")+ylab("")+ggtitle({.y})+
         scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
         theme(plot.title  = element_text(size = 9,face="italic"),
               panel.border = element_blank(),
               axis.title = element_text(size=9,face="plain"))+ 
         scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals))

# add y label
dif5_its[[1]] <-  dif5_its[[1]]+ylab("URA (%)") 

# add gap for outlier
dif5_its[[1]] <- gg.gap( dif5_its[[1]], ylim = c(4,100), segments = c(6,90), tick_width = c(2,2),rel_heights = c(1,0,2.5))
dif5_its[[1]]<- ggdraw(dif5_its[[1]])+draw_label("//", x=0.11,y=0.305, angle=305)
dif5_its <- plot_grid(plotlist = dif5_its , nrow=1)

# extract coefficients
tss_coefficients_its <- its_difs[[5]]$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Must_TSS")) %>% 
  mutate(rowname=its_difs[[5]]$significant_taxa)

# combine plots
fung_difabund <- dif1_its/dif2_its/dif4_its/dif5_its+plot_annotation(tag_levels = "A")


#WriteXLS(x=c("vintage_coefficients_its", "precipitation_coefficients_its",  "ph_coefficients_its","tss_coefficients_its"), ExcelFileName = "Table_S6.xls")
ggsave("Figures_and_tables_check/Figure_7_fungi_differential_abundance.pdf", fung_difabund, width = 8.5,height = 11, units = "in"  )

