#library(vegan)
#library(microbiome)
#library(microbiomeSeq)
#library(RColorBrewer)
#library(phyloseq)
#library(shades)
#library(WriteXLS)
#source("setup.R")

# alpha diversity function
plot_alpha_div <- function(physeq){
  set.seed(2)
  rdf <- physeq %>% rarefy_even_depth() %>% otu_table() %>% as.data.frame 
hills <-  rdf %>%
    renyi(hill=TRUE, scales = c(0,1,2)) %>% rownames_to_column("Shipment") %>% rename("Richness"=`0`, "Exp. Shannon"=`1`, "Inverse Simpson"=`2`) %>% left_join(merge_samdat, by="Shipment") %>% 
    pivot_longer(`Richness`:`Inverse Simpson`, names_to = "Diversity")%>% mutate(Diversity=factor(Diversity, levels=c("Richness", "Exp. Shannon", "Inverse Simpson"), ordered = TRUE)) %>% group_by(Year,Diversity)%>%
    mutate(meanval=mean(value))
  
  figure <- hills %>% ggplot(aes(group=AVA, fill=AVA, shape=AVA,y=value,x=AVA))+
    geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
    geom_jitter(aes(),size=3,width=0.2,height=0, show.legend = FALSE)+ 
    facet_grid(rows = vars(Diversity), cols=vars(Year), scales = "free_y", switch = "y")+
    geom_hline(mapping = aes(yintercept=meanval),linetype="dotted", color="black", size=1, alpha=0.1)+
    scale_fill_manual(values=mypal)+scale_shape_manual(values=shapevals)+
      theme(strip.placement = "outside", strip.background = element_rect(fill="white", color="white"))+
    theme(axis.text.x = element_text(angle=45, hjust=1))+xlab("")+ylab("")
  return(list(figure,hills))
  
}



## Supplemental figure 1: Bacterial alpha diversity

bact_alpha <- plot_alpha_div(must_16s_merge)

bact_alpha_test_df <- bact_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment ) %>% left_join(merge_samdat)

bact_alpha_long <- read.csv("bact_alpha_long.csv")

bacterial_alpha_diversity <- bact_alpha_long %>% 
  mutate(Diversity=factor(Diversity, levels=c("Richness", "Exp. Shannon", "Inverse Simpson"), ordered = TRUE)) %>%
  ggplot(aes(group=AVA, fill=AVA, shape=AVA,y=value,x=AVA))+
  geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
  geom_jitter(aes(),size=3,width=0.2,height=0, show.legend = FALSE)+ 
  geom_text(aes(x=AVA, y=mean_sep_y, label=mean_sep_val))+
  facet_grid(rows = vars(Diversity), cols=vars(Year), scales = "free_y", switch = "y")+
  geom_hline(mapping = aes(yintercept=meanval),linetype="dotted", color="black", size=1, alpha=0.4)+
  scale_fill_manual(values=mypal)+scale_shape_manual(values=shapevals)+
  theme(strip.placement = "outside", strip.background = element_rect(fill="white", color="white"))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+xlab("")+ylab("")


# commented out because assumptions and descriptive stats uses this script to generate alpha div values and don't want it to print out extraneous figs.
#ggsave( "Figures_and_tables_check/bacterial_alpha_diversity.pdf",bacterial_alpha_diversity,width = 6.5,height = 6.5, units = "in" )


## Supplemental figure 2: Fungal alpha diversity

fung_alpha <- plot_alpha_div(must_its_merge)
Supplemental_figure_2_fungal_alpha_diversity <- fung_alpha[[1]]

fungt_alpha_test_df <- fung_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment ) %>% left_join(merge_samdat)

Supplemental_figure_2_fungal_alpha_diversity

fung_alpha_long <- read.csv("Files/fung_alpha_long.csv")

fungal_alpha_diversity <- fung_alpha_long %>% 
 mutate(Diversity=factor(Diversity, levels=c("Richness", "Exp. Shannon", "Inverse Simpson"), ordered = TRUE)) %>%
  ggplot(aes(group=AVA, fill=AVA, shape=AVA,y=value,x=AVA))+
  geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
  geom_jitter(aes(),size=3,width=0.2,height=0, show.legend = FALSE)+ 
  geom_text(aes(x=AVA, y=mean_sep_y, label=mean_sep_val))+
  facet_grid(rows = vars(Diversity), cols=vars(Year), scales = "free_y", switch = "y")+
  geom_hline(mapping = aes(yintercept=meanval),linetype="dotted", color="black", size=1, alpha=0.4)+
  scale_fill_manual(values=mypal)+scale_shape_manual(values=shapevals)+
  theme(strip.placement = "outside", strip.background = element_rect(fill="white", color="white"))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+xlab("")+ylab("")


# commented out because assumptions and descriptive stats uses this script to generate alpha div values and don't want it to print out extraneous figs.
#ggsave("Figures_and_tables_check/fungal_alpha_diversity.pdf",fungal_alpha_diversity, width=6.5,height=6.5, units="in")

#let's also save the tables, just for reference
#fung_alpha_df <- fung_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment )
#bact_alpha_df <- bact_alpha[[2]] %>% pivot_wider(names_from=Diversity, id_cols = Shipment )
#WriteXLS(ExcelFileName = "Figures_and_tables_check/alpha_diversities.xls", x=c("fung_alpha_df","bact_alpha_df"))



# rarefaction curves

must_16s_merge_inext <- must_16s_merge
otu_table(must_16s_merge_inext) <- otu_table(must_16s_merge_inext) %>% t
#inext0_bact <- phyloseq_inext(must_16s_merge_inext)
inext0_bact <- readRDS("bact_rarecurve_rich.rds")+ylab("Richness")
inext0_bact<- remove_geom(inext0_bact, "label")

#inext1_bact <- phyloseq_inext(must_16s_merge_inext, Q=1, multithread = TRUE)
inext1_bact <- readRDS("bact_rarecurve_expshan.rds")+ylab("Exp. Shannon")
inext1_bact<- remove_geom(inext1_bact, "label")

#inext2_bact <- phyloseq_inext(must_16s_merge_inext, Q=2, multithread = TRUE)
inext2_bact <- readRDS("bact_rarecurve_invsimp.rds")+ylab("Inv. Simpson")
inext2_bact<- remove_geom(inext2_bact, "label")


must_its_merge_inext <- must_its_merge
otu_table(must_its_merge_inext) <- otu_table(must_its_merge_inext) %>% t

#inext0_fung<- phyloseq_inext(must_its_merge_inext)
inext0_fung <- readRDS("fung_rarecurve_rich.rds") + ylab("Richness")
inext0_fung<- remove_geom(inext0_fung, "label")

#inext1_fung <- phyloseq_inext(must_its_merge_inext, Q=1, multithread = TRUE)
inext1_fung <- readRDS("fung_rarecurve_expshan.rds")+ ylab("Exp. Shannon")
inext1_fung<- remove_geom(inext1_fung, "label")

#inext2_fung <- phyloseq_inext(must_its_merge_inext, Q=2, multithread = TRUE)
inext2_fung <- readRDS("fung_rarecurve_invsimp.rds")+ylab("Inv. Simpson")
inext2_fung<- remove_geom(inext2_fung, "label")


bact_rarecurve <- inext0_bact/inext1_bact/inext2_bact+plot_annotation(tag_levels = "A")
#ggsave(plot = bact_rarecurve,filename  =  "Supplemental_figure_1_rarefaction_bacteria.pdf", device = "pdf", width=6.5, height=9, units="in")

fung_rarecurve <- inext0_fung/inext1_fung/inext2_fung+plot_annotation(tag_levels = "A")
#ggsave(plot=fung_rarecurve, "Supplemental_figure_2_rarefaction_fungi.pdf", device="pdf",  width=6.5, height=9, units="in")
