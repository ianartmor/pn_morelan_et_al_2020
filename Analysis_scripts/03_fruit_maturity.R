#fruit maturity figure

#library(RColorBrewer)
#library(shades)
#library(tidyverse) 
#source("setup.R")
fm_df <- merge_samdat %>%
  pivot_longer(cols = c(Must_TSS , Must_pH, Must_TA), names_to = "Chem_var") %>% 
  group_by(Year,Chem_var)%>%
  mutate(meanval=mean(value))

fm_df$facets = factor(fm_df$Chem_var, labels = c("TA~(g~tartaric~acid~L^{-1})", "pH", "TSS~('°Brix')"))

Figure_3_fruit_maturity <- fm_df %>%  group_by(AVA,Year) %>% 
  ggplot(aes(x=AVA,y=value,fill=AVA))+
  geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
  facet_grid(rows=vars(facets), cols=vars(Year), scales = "free_y", labeller = "label_parsed", switch = "y")+
  geom_hline(mapping = aes(yintercept=meanval),linetype="dotted", color="black", size=1)+
  geom_jitter(size=2, aes(shape=AVA), width=0.2,height=0, show.legend = FALSE)+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")+
  ylab("")+xlab("")+ scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals)+ 
  theme(strip.placement = "outside", strip.background = element_rect(fill="white", color="white"))

Figure_3_fruit_maturity
#ggsave("Figures_and_tables_check/Figure_3_fruit_maturity.pdf",Figure_3_fruit_maturity, width = 6.5,height = 9, units = "in" )
