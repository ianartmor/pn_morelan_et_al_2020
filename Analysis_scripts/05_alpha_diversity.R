# alpha diversity

plot_alpha_div <- function(physeq){
  
  hills <- physeq %>% rarefy_even_depth() %>% otu_table() %>% as.data.frame %>%
    renyi(hill=TRUE, scales = c(0,1,2)) %>% rownames_to_column("Shipment") %>% rename("Richness"=`0`, "Exp. Shannon"=`1`, "Inverse Simpson"=`2`) %>% left_join(merge_samdat, by="Shipment") %>% 
    pivot_longer(`Richness`:`Inverse Simpson`, names_to = "Diversity")%>% mutate(Diversity=factor(Diversity, levels=c("Richness", "Exp. Shannon", "Inverse Simpson"), ordered = TRUE)) %>% group_by(Year,Diversity)%>%
    mutate(meanval=mean(value))
  
  figure <- hills %>% ggplot(aes(group=AVA, fill=AVA, shape=AVA,y=value,x=AVA))+
    geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
    geom_jitter(aes(),size=3,width=0.2,height=0, show.legend = FALSE)+ 
    facet_grid(rows = vars(Diversity), cols=vars(Year), scales = "free_y")+
    geom_hline(mapping = aes(yintercept=meanval),linetype="dotted", color="black", size=1)+
    scale_fill_manual(values=mypal)+scale_shape_manual(values=shapevals)+
    theme(axis.text.x = element_text(angle=45, hjust=1))+xlab("")+ylab("Diversity")
  return(list(figure,hills))
  
}



## Supplemental figure 1: Bacterial alpha diversity

bact_alpha <- plot_alpha_div(must_16s_merge)

bact_alpha_test_df <- bact_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment ) %>% left_join(merge_samdat)

wilcox.test(
  x=bact_alpha_test_df %>% filter(Year ==2016) %>% select(Richness)%>% unlist %>% as.numeric(), 
  y=bact_alpha_test_df %>% filter(Year ==2017) %>% select(Richness) %>% unlist %>% as.numeric())

wilcox.test(
  x=bact_alpha_test_df %>% filter(Year ==2016) %>% select(`Exp. Shannon`)%>% unlist %>% as.numeric(), 
  y=bact_alpha_test_df %>% filter(Year ==2017) %>% select(`Exp. Shannon`) %>% unlist %>% as.numeric())

wilcox.test(
  x=bact_alpha_test_df %>% filter(Year ==2016) %>% select(`Inverse Simpson`)%>% unlist %>% as.numeric(), 
  y=bact_alpha_test_df %>% filter(Year ==2017) %>% select(`Inverse Simpson`) %>% unlist %>% as.numeric())

kruskal.test(x=bact_alpha_test_df$Richness, g=bact_alpha_test_df$AVA)


kruskal.test(x=bact_alpha_test_df$`Exp. Shannon`, g=bact_alpha_test_df$AVA)
kruskal.test(x=bact_alpha_test_df$`Inverse Simpson`, g=bact_alpha_test_df$AVA)

Supplemental_figure_1_bacterial_alpha_diversity <- bact_alpha[[1]]

Supplemental_figure_1_bacterial_alpha_diversity
ggsave( "Figures_and_tables_check/Supplemental_figure_1_bacterial_alpha_diversity.pdf",Supplemental_figure_1_bacterial_alpha_diversity)

fung_alpha <- plot_alpha_div(must_its_merge)
Supplemental_figure_2_fungal_alpha_diversity <- fung_alpha[[1]]

fungt_alpha_test_df <- fung_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment ) %>% left_join(merge_samdat)

wilcox.test(
  x=fungt_alpha_test_df %>% filter(Year ==2016) %>% select(Richness)%>% unlist %>% as.numeric(), 
  y=fungt_alpha_test_df %>% filter(Year ==2017) %>% select(Richness) %>% unlist %>% as.numeric())

wilcox.test(
  x=fungt_alpha_test_df %>% filter(Year ==2016) %>% select(`Exp. Shannon`)%>% unlist %>% as.numeric(), 
  y=fungt_alpha_test_df %>% filter(Year ==2017) %>% select(`Exp. Shannon`) %>% unlist %>% as.numeric())

wilcox.test(
  x=fungt_alpha_test_df %>% filter(Year ==2016) %>% select(`Inverse Simpson`)%>% unlist %>% as.numeric(), 
  y=fungt_alpha_test_df %>% filter(Year ==2017) %>% select(`Inverse Simpson`) %>% unlist %>% as.numeric())

kruskal.test(x=fungt_alpha_test_df$Richness, g=fungt_alpha_test_df$AVA)
kruskal.test(x=fungt_alpha_test_df$`Exp. Shannon`, g=fungt_alpha_test_df$AVA)
kruskal.test(x=fungt_alpha_test_df$`Inverse Simpson`, g=fungt_alpha_test_df$AVA)

Supplemental_figure_2_fungal_alpha_diversity
ggsave("Figures_and_tables_check/Supplemental_figure_2_fungal_alpha_diversity.pdf",Supplemental_figure_2_fungal_alpha_diversity)

#let's also save the tables, just for reference
fung_alpha_df <- fung_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment )
bact_alpha_df <- bact_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment )
WriteXLS(ExcelFileName = "Figures_and_tables_check/alpha_diversities.xls", x=c("fung_alpha_df","bact_alpha_df"))
