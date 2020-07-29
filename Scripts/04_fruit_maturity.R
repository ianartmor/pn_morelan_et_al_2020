Figure_3_fruit_maturity <- merge_samdat %>%
  rename("TSS (°Brix)"=Must_TSS, "pH"=Must_pH, "TA (g/L)"=Must_TA) %>%
  pivot_longer(cols = c(`TSS (°Brix)` , `pH`, `TA (g/L)`), names_to = "Chem_var") %>% 
  group_by(Year,Chem_var)%>%
  mutate(meanval=mean(value))%>%
  group_by(AVA,Year) %>% 
  ggplot(aes(x=AVA,y=value,fill=AVA))+
  geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
  facet_grid(rows=vars(Chem_var), cols=vars(Year), scales = "free_y")+
  geom_hline(mapping = aes(yintercept=meanval),linetype="dotted", color="black", size=1)+
  geom_jitter(size=3, aes(shape=AVA), width=0.2,height=0, show.legend = FALSE)+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")+
  ylab("")+xlab("")+ scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals)

Figure_3_fruit_maturity
ggsave("Figures_and_tables_check/Figure_3_fruit_maturity.pdf",Figure_3_fruit_maturity)