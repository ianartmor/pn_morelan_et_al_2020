#library(dada2)
#library(microbiomeSeq)
#library(microbiome)
#library(corncob)
#library(patchwork)
#library(cowplot)
#library(tidyverse)
#library(phyloseq)
#source("Analysis_scripts/00_setup.R")
#source("Analysis_scripts/07a_da_functions.R")



# differential abundance using corncob

set.seed(5)
  

# extract tax table
taxtab <- must_its_merge %>% tax_table() %>% data.frame

#replace "g__"
genus_names <- taxtab$Genus %>% as.character() %>% gsub("g__", "", .)

#update tax table
taxtab$Genus <- genus_names
tax_table(must_its_merge) <- tax_table(taxtab %>% as.matrix)

its_merge_gen <- taxa_level(must_its_merge, "Genus")
toptax <- prune_taxa(top_taxa(its_merge_gen, 8),its_merge_gen)

#outlier_detect(toptax)

#no_outliers_its <- toptax %>% 
#  subset_samples(Shipment!="AS1_2017") %>% 
#  subset_samples(Shipment!="RRV1_2016") %>% 
#  subset_samples(Shipment!="RRV1_2017") %>% 
#  subset_samples(Shipment != "SNC1_2016") %>% 
#  subset_samples(Shipment!="CRN1_2017") 

### VINTAGE

vintage_res_its <- vintage_dif(toptax)

#toptax
sig_vintage_its <- toptax  %>%  transform("compositional") %>%  prune_taxa(taxa = vintage_res_its$significant_taxa,x=.)


#extract coefficients
vintage_coefficients_its <- vintage_res_its$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Year")) %>% 
  mutate(rowname=vintage_res_its$significant_taxa)

#toptax and plot
vintage_plot <- sig_vintage_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
  purrr::map2(.x = ., .y=names(.), .f = ~ ggplot(.,aes(x=as.factor(Year), y=Abundance*100, group=Year, fill=AVA, shape=as.factor(Year)))+
                geom_point(show.legend = FALSE)+geom_boxplot(alpha=0.2, show.legend = FALSE)+
                xlab("")+ylab("")+ggtitle({.y})+
                theme(plot.title  = element_text(size = 9, face="italic"),
                      axis.title = element_text(size=9, face="plain"), 
                      panel.border = element_blank(),
                     # legend.title = element_text(size=9),
                      #legend.text = element_text(size=7), legend.spacing.y = unit(1, "mm"), legend.box = "horizontal"
                      )+ 
                scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals, name="Vintage")+
                guides(fill = guide_legend(override.aes=list(shape=21)))+ 
                scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
                scale_fill_manual(values = mypal))

vintage_plot[[1]] <- vintage_plot[[1]] + ylab("URA (%)")
vintage_plots <- plot_grid(plotlist = vintage_plot, nrow=1)

### PRECIPITATION
precip_res_its <- must_precip_growing_dif(toptax)

# toptax
must_precip_sig_its <- toptax%>%
  transform("compositional") %>%  prune_taxa(precip_res_its$significant_taxa,.)

precip_coefficients_its <- precip_res_its$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Precipitation.growing")) %>% 
  mutate(rowname=precip_res_its$significant_taxa)

#split and plot
precip_plot <- must_precip_sig_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
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
precip_plot[[1]] <-  precip_plot[[1]]+ylab("URA (%)") 

# add gaps for outliers
precip_plot[[2]] <-  gg.gap( precip_plot[[2]], ylim = c(0,94), segments = c(10,92), tick_width = c(2,2), rel_heights = c(3,0,1))
precip_plot[[2]]<- ggdraw(precip_plot[[2]])+draw_label("//", x=0.31,y=0.755, angle=305)

precip_plot[[3]] <- gg.gap( precip_plot[[3]], ylim = c(4,100), segments = c(6,90), tick_width = c(2,2),rel_heights = c(1,0,2.5))
precip_plot[[3]]<- ggdraw(precip_plot[[3]])+draw_label("//", x=0.335,y=0.295, angle=305)

precip_plots <- plot_grid(plotlist = precip_plot, nrow=1)

### TA
# no significant taxa
#ta_res <- must_ta_dif(toptax)

#### PH

ph_res_its <- must_pH_dif(toptax)

#toptax
must_ph_sig_its <- toptax%>% 
  transform("compositional") %>%  prune_taxa(ph_res_its$significant_taxa,.)

#extract coefficients
ph_coefficients_its <- ph_res_its$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,  .f= ~filter(.,rowname=="mu.Must_pH")) %>% 
  mutate(rowname=ph_res_its$significant_taxa)

# split and plot
ph_plot <- must_ph_sig_its %>%psmelt %>% filter(!grepl("f__", OTU)) %>% split(.$OTU) %>% 
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
ph_plot[[1]] <-  ph_plot[[1]]+ylab("URA (%)") 

#add gap for outlier
ph_plot[[1]] <- gg.gap( ph_plot[[1]], ylim = c(4,100), segments = c(6,90), tick_width = c(2,2),rel_heights = c(1,0,2.5))
ph_plot[[1]]<- ggdraw(ph_plot[[1]])+draw_label("//", x=0.43,y=0.295, angle=305)

ph_plots <- plot_grid(plotlist=ph_plot, nrow=1)




### TSS

tss_res_its <- must_tss_dif(toptax)
#toptax
must_tss_sig_its <- toptax %>%
  transform("compositional") %>%  prune_taxa(tss_res_its$significant_taxa,.)

#split and plot
tss_plot <- must_tss_sig_its %>%psmelt %>% mutate(Vintage=as.factor(Year)) %>%filter(!grepl("f__", OTU)) %>%split(.$OTU) %>% 
  map2(.x=.,.y=names(.),.f=~ggplot(.,aes(x=Must_TSS, y=Abundance*100, fill=AVA, shape=Vintage))+
         geom_point()+ 
         # geom_smooth(method="lm", se=FALSE)+
         xlab("TSS (°Brix)")+ylab("")+ggtitle({.y})+
         scale_y_continuous(labels = scales::number_format(accuracy = 0.01))+
         theme(plot.title  = element_text(size = 9,face="italic"),
               panel.border = element_blank(),
               axis.title = element_text(size=9,face="plain"))+ 
         scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals))

# add y label
tss_plot[[1]] <-  tss_plot[[1]]+ylab("URA (%)") 

# add gap for outlier
tss_plot[[1]] <- gg.gap(tss_plot[[1]], ylim = c(4,100), segments = c(6,90), tick_width = c(2,2),rel_heights = c(1,0,2.5))
tss_plot[[1]]<- ggdraw(tss_plot[[1]])+draw_label("//", x=0.11,y=0.305, angle=305)

tss_plots <- plot_grid(plotlist = tss_plot , nrow=1)

# extract coefficients
tss_coefficients_its <- tss_res_its$significant_models %>% purrr::map(coef) %>% map(data.frame) %>%
  map(rownames_to_column) %>% purrr::map_df(.,.f=~filter(.,rowname=="mu.Must_TSS")) %>% 
  mutate(rowname=tss_res_its$significant_taxa)

# combine plots
fung_difabund <- vintage_plots/precip_plots/ph_plots/tss_plots+plot_annotation(tag_levels = "A")


WriteXLS(x=c("vintage_coefficients_its", "precip_coefficients_its",  "ph_coefficients_its","tss_coefficients_its"), ExcelFileName = "Table_S6_corncob_fung.xls")
ggsave("Figures_and_tables_check/Figure_7_fungi_differential_abundance.pdf", fung_difabund, width = 8.5,height = 11, units = "in"  )

