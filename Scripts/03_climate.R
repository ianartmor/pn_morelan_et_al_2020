Figure_2_climate <- merge_samdat %>%
  rename( "Precipitation (mm)"=Precipitation.growing, "Growing Degree Days"=Growing.Degree.Days) %>%
  pivot_longer(cols = c(`Precipitation (mm)`, `Growing Degree Days`), names_to = "Climate_var") %>% 
  group_by(Year,Climate_var)%>%
  mutate(meanval=mean(value))%>%
  group_by(AVA,Year) %>% 
  ggplot(aes(x=AVA,y=value,fill=AVA))+
  geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
  facet_grid(rows=vars(Climate_var), cols=vars(Year), scales = "free_y")+
  geom_hline(mapping = aes(yintercept=meanval),linetype="dotted")+
  geom_jitter(aes(shape=AVA),size=3,width=0.2,height=0,show.legend = FALSE)+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")+
  ylab("")+xlab("")+
  scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals)

kruskal.test(x=merge_samdat %>% filter(Year==2016) %>% select(Growing.Degree.Days) %>% unlist, g=merge_samdat %>% filter(Year==2016) %>% select(AVA) %>% unlist)

kruskal.test(x=merge_samdat %>% filter(Year==2017) %>% select(Growing.Degree.Days) %>% unlist, g=merge_samdat %>% filter(Year==2017) %>% select(AVA) %>% unlist)

kruskal.test(x=merge_samdat %>% filter(Year==2016) %>% select(Precipitation.growing) %>% unlist, g=merge_samdat %>% filter(Year==2016) %>% select(AVA) %>% unlist)

dunn.test(x=merge_samdat %>% filter(Year==2016) %>% select(Precipitation.growing) %>% unlist, g=merge_samdat %>% filter(Year==2016) %>% select(AVA) %>% unlist)


kruskal.test(x=merge_samdat %>% filter(Year==2017) %>% select(Precipitation.growing) %>% unlist, g=merge_samdat %>% filter(Year==2017) %>% select(AVA) %>% unlist)


dunn.test(x=merge_samdat %>% filter(Year==2017) %>% select(Precipitation.growing) %>% unlist, g=merge_samdat %>% filter(Year==2017) %>% select(AVA) %>% unlist)


wilcox.test(
  x=merge_samdat %>% filter(Year ==2016) %>% select(Growing.Degree.Days)%>% unlist %>% as.numeric(), 
  y=merge_samdat %>% filter(Year ==2017) %>% select(Growing.Degree.Days) %>% unlist %>% as.numeric())

wilcox.test(
  x=merge_samdat %>% filter(Year ==2016) %>% select(Precipitation.growing)%>% unlist %>% as.numeric(), 
  y=merge_samdat %>% filter(Year ==2017) %>% select(Precipitation.growing) %>% unlist %>% as.numeric())

Figure_2_climate
ggsave("Figures_and_tables_check/Figure_2_climate.pdf" ,Figure_2_climate)