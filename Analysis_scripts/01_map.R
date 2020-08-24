
# download maps
# google map API key required
ca_or <- get_map(location = "Sacramento", zoom = 5, maptype = "terrain-background")
norcal_zoom <- get_map(location = "Sonoma", zoom = 8, maptype = "terrain-background")

# Calculate average Precip, GDD, and combine with relevant sample data for mapping
map_df <-merge_samdat %>% 
  group_by(Vineyard) %>% 
  summarise(Precipitation.growing = mean(Precipitation.growing), Growing.Degree.Days= mean(Growing.Degree.Days)) %>% 
  left_join(merge_samdat %>% 
              select(AVA, Vineyard, Latitude, Longitude), by= "Vineyard") %>% 
  unique()

#plot zoomed in section
norcal_plot <- ggmap(norcal_zoom, extent="panel")+
  # Color is precipitation, size is GDD
  geom_point(aes(x=Longitude,y=Latitude,color= Precipitation.growing, size=Growing.Degree.Days), data = map_df)+
  # manually change color gradient
  scale_color_gradient2(name="Precipitation (mm)", low = "red", mid="gray40", high="blue", midpoint = mean(map_df$Precipitation.growing))+
  # repel labels from points for readability
  geom_label_repel(data=map_df, aes(label=Vineyard,x=Longitude,y=Latitude, fill=AVA), show.legend = FALSE, xlim=c(-122,-121))+labs(x="", y="")+
  scale_fill_manual(values = mypal)+ 
  scale_size_continuous(name="Growing Degree Days")+
  coord_equal()+
  
  ggsn::scalebar(y.min = 37.25, y.max = 38,
                 x.min = -124, x.max = -123,  dist = 25, dist_unit = "km",
                 transform  = TRUE, model = 'WGS84', height=0.1,st.dist=0.15, st.size = 2)+
  
  ggsn::north(y.min = 36.9, y.max = 37.75,
              x.min = -124, x.max = -123.75, scale = 1 ,symbol = 3
  )


# extract values of border of above plot so a gray box can be plotted on larger plot
box <-data.frame(vals= c(norcal_plot$scales$scales[[1]]$limits, norcal_plot$scales$scales[[2]]$limits))

# plot zoomed out map
ca_or_plot <- ggmap(ca_or, extent = "panel")+
  # plot the gray box
  geom_rect(inherit.aes = FALSE, mapping = aes(xmin=box$vals[1], xmax=box$vals[2], ymin=box$vals[3],ymax=box$vals[4]), alpha = 0.1)+
  # same point mapping as first plot
  geom_point(aes(x=Longitude,y=Latitude,color= Precipitation.growing, size=Growing.Degree.Days), data = map_df)+
  scale_color_gradient2(name= "Precipitation (mm)",low = "red", mid="gray40", high="blue", midpoint = mean(map_df$Precipitation.growing))+
  geom_label_repel(data=map_df, aes(label=Vineyard,x=Longitude,y=Latitude, fill=AVA), show.legend = FALSE, xlim=c(-127,-125))+
  scale_x_continuous(limits = c(-127,-117), name = "")+
  scale_y_continuous(limits=c(32,46), name="")+
  theme(legend.position = "none")+
  scale_fill_manual(values = mypal)+coord_equal()+
  
  ggsn::scalebar(y.min = 32.5, y.max = 33.5,
                 x.min = -120, x.max = -117,  dist = 100, dist_unit = "km",
                 transform  = TRUE, model = 'WGS84', height=0.25,st.dist=0.5, st.size = 2)+
  
  ggsn::north(y.min = 32.5, y.max = 33.5,
              x.min = -120.5, x.max = -118.5, scale = 1.5,symbol = 3
  )


#combine and arrange plots using patchwork
Figure_1_map <- ca_or_plot|(norcal_plot/guide_area()+plot_layout(guides="collect", heights=c(4,1)) & 
                              theme( legend.box="horizontal", 
                                     plot.margin = margin(0,0,0,0,"cm"))) + 
  plot_layout(widths=c(1,2))& 
  theme(plot.margin = margin(0,0,0,0,"cm"))  

Figure_1_map <- Figure_1_map+plot_annotation(tag_levels = "A")& 
  theme(plot.tag.position = c(0.05, 0.95))

Figure_1_map
#save as PDF
ggsave("Figures_and_tables_check/Figure_1_map.pdf",Figure_1_map,  width = 8, height = 9, units = "in")