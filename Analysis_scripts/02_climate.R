#library(tidyverse)
#library(RColorBrewer)
#library(shades)
#source("setup.R")

#climate_long <- merge_samdat %>%
 # rename( "Precipitation (mm)"=Precipitation.growing, "Growing Degree Days"=Growing.Degree.Days) %>%
  #pivot_longer(cols = c(`Precipitation (mm)`, `Growing Degree Days`), names_to = "Climate_var") 
climate_long <- read.csv("Data/climate_long.csv")

Figure_2_climate <- climate_long %>% 
  group_by(Year,Climate_var)%>%
  mutate(meanval=mean(value))%>%
  group_by(AVA,Year) %>% 
  ggplot(aes(x=AVA,y=value,fill=AVA))+
  geom_boxplot(alpha=0.1, color="gray70", show.legend = FALSE)+
  geom_text(aes(x=AVA, y=mean_sep_y, label=mean_sep_val))+
  facet_grid(rows=vars(Climate_var), cols=vars(Year), scales = "free_y", switch = "y")+
  geom_hline(mapping = aes(yintercept=meanval),linetype="dotted")+
  geom_jitter(aes(shape=AVA),size=2,width=0.2,height=0,show.legend = FALSE)+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")+
  ylab("")+xlab("")+
  scale_fill_manual(values = mypal)+scale_shape_manual(values=shapevals)+ 
  theme(strip.placement = "outside", strip.background = element_rect(fill="white", color="white"))

Figure_2_climate
#ggsave("Figures_and_tables_check/Figure_2_climate.pdf" ,Figure_2_climate,width = 6.5,height = 6, units = "in" )
