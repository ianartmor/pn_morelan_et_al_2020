#beta diversity

source("setup.R")

#check if replicates are different. 
#must_16s %>% rarefy_even_depth(sample.size = 5000)%>% ordinate(method="NMDS") %>% plot_ordination(physeq=must_16s, ordination = .)+facet_wrap(~Shipment) 

#beta diversity function
plot_beta_div <- function(physeq_merge){
  
  merge_samdat <- sample_data(physeq_merge) %>%data.frame 

    set.seed(1)
  results_perm <- adonis(formula=otu_table(physeq_merge)~Year+AVA, method = "bray",
                         data=merge_samdat, permutations = 9999)
  
  #physeq.ord <- capscale(formula = transform(otu_table(top20),"hellinger")~1, 
                      #   distance = "euclidean")
  physeq_full.ord <- metaMDS(otu_table(physeq_merge), k=3,
                             distance = "bray")
  
  to_plot <- physeq_full.ord$points %>% as.data.frame()%>% rownames_to_column("Shipment") %>% 
    left_join(merge_samdat, by="Shipment")
  set.seed(1)
  
  results_envfit <- envfit(physeq_full.ord, merge_samdat, permutations = 9999 )
  
  bdiv_ava <- gg_ordiplot(physeq_full.ord, groups = merge_samdat$AVA, 
                          spiders = FALSE, hull = FALSE, ellipse = FALSE, plot=FALSE)
  
  bdiv_ava[[6]]$layers[[1]] <- NULL
  bdiv_ava$df_ord$Year <- bdiv_ava$df_ord %>% rownames %>% str_sub(-4) %>% as.numeric()
  bdiv_plot <- bdiv_ava[[6]]+  
    scale_color_manual(values=mypal)+
    geom_point(data=bdiv_ava$df_ord,aes(shape=Group,x=x,y=y,fill=Group, size=as.factor(Year)))+
    scale_shape_manual(values=shapevals, name="AVA")+
    scale_fill_manual(values=mypal, name="AVA")+
    scale_size_discrete(range=c(2,4), name="Vintage")+
     scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    xlab("NMDS1")+ylab("NMDS2")+coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))
  
  sample_plot <- bdiv_plot+coord_fixed()
  
  bdiv_env <- ggordiplots::gg_envfit(physeq_full.ord, groups = merge_samdat$AVA, env=merge_samdat %>% select(-Longitude, "pH"=Must_pH, "TSS"=Must_TSS,"TA"=Must_TA), arrow.col = "black", plot=FALSE)
  df_arrows <- bdiv_env$df_arrows
  
  vector_plot1 <-gg_ordisurf(ord = physeq_full.ord, env.var = merge_samdat$Precipitation.dormant, plot = FALSE)
  vector_plot1 <- vector_plot1$plot+ 
    xlab("NMDS1")+ylab("NMDS2")+
    coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    scale_color_continuous(name="Precipitation (mm)")+
    xlab("NMDS1")+ylab("NMDS2")+coord_fixed()
  
  
  figure<- sample_plot / vector_plot1 # +plot_annotation(tag_levels = "A")
  
  return (
    list(
      figure,
      results_perm,
      results_envfit
    )
  )
}

# mantel function
physeq_mantel <- function(physeq){
  set.seed(1)
  otutab <- otu_table(physeq) %>% transform(transform="hellinger") %>% as.matrix()
  
  xy <- physeq %>% sample_data%>% as_tibble()%>% select("y"= Latitude, "x"=Longitude) 
  
  xydist <- geodist(xy)
  bcdist <- vegdist(otutab, "euclidean")
  mantel(xydist, bcdist,permutations = 1000, method = "spearman")
}

# run beta div function
beta_div_16s <- plot_beta_div(must_16s_merge)

#extract figure
Figure_4_bacterial_beta_diversity <- beta_div_16s[[1]]

# add annotation with patchwork
Figure_4_bacterial_beta_diversity <- Figure_4_bacterial_beta_diversity+plot_annotation(tag_levels = "A", title = "Bacterial Communities")&
  theme(plot.title = element_text(hjust = 0.3))
ggsave("Figures_and_tables_check/Figure_4_bacterial_beta_diversity.pdf", Figure_4_bacterial_beta_diversity, width = 8,height = 8,units = "in" )

#mantel tests
physeq_mantel(must_16s_2016)
physeq_mantel(must_16s_2017)

#run beta div function
beta_div_its <- plot_beta_div(must_its_merge)

#extract figure
Figure_5_fungal_beta_diversity <- beta_div_its[[1]]

#add annotation with patchwork
Figure_5_fungal_beta_diversity <- Figure_5_fungal_beta_diversity+plot_annotation(tag_levels = "A", title = "Fungal Communities")&
  theme(plot.title = element_text(hjust = 0.3))
ggsave("Figures_and_tables_check/Figure_5_fungal_beta_diversity.pdf", Figure_5_fungal_beta_diversity, width = 8,height = 8,units = "in"  )

#mantel tests
physeq_mantel(must_its_2016)
physeq_mantel(must_its_2017)




