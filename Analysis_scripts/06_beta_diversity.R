#beta diversity


plot_beta_div <- function(physeq_merge){
  
  merge_samdat <- sample_data(physeq_merge) %>%data.frame 
  top20 <- physeq_merge %>% prune_taxa(top_taxa(.,20),.)
  set.seed((1))
  results_perm <- adonis(formula =transform(otu_table(physeq_merge),"hellinger")~Year+AVA, data=merge_samdat, permutations = 9999)
  physeq.ord <- capscale(formula = transform(otu_table(top20),"hellinger")~1, distance = "euclidean")
  physeq_full.ord <- metaMDS(transform(otu_table(physeq_merge),"hellinger"), k=3)
  
  to_plot <- physeq.ord$points %>% as.data.frame()%>% rownames_to_column("Shipment") %>% left_join(merge_samdat, by="Shipment")
  set.seed(1)
  
  results_envfit <- envfit(physeq.ord, merge_samdat )
  
  bdiv_ava <- gg_ordiplot(physeq.ord, groups = merge_samdat$AVA, spiders = TRUE, hull = FALSE, ellipse = FALSE, plot=FALSE)
  bdiv_ava[[6]]$layers[[1]] <- NULL
  bdiv_ava$df_ord$Year <- bdiv_ava$df_ord %>% rownames %>% str_sub(-4) %>% as.numeric()
  bdiv_plot <- bdiv_ava[[6]]+  
    scale_color_manual(values=mypal)+
    geom_point(data=bdiv_ava$df_ord,aes(shape=Group,x=x,y=y,fill=Group, size=as.factor(Year)))+
    scale_shape_manual(values=shapevals, name="AVA")+
    scale_fill_manual(values=mypal, name="AVA")+
    scale_size_discrete(range=c(2,4), name="Year")+
    xlab("NMDS1")+ylab("NMDS2")+coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))
  
  sample_plot <- bdiv_plot+coord_fixed()
  
  bdiv_env <- ggordiplots::gg_envfit(physeq_full.ord, groups = merge_samdat$AVA, env=merge_samdat %>% select(-Longitude, "pH"=Must_pH, "TSS"=Must_TSS,"TA"=Must_TA), arrow.col = "black", plot=FALSE)
  df_arrows <- bdiv_env$df_arrows
  
  vector_plot1 <-gg_ordisurf(ord = physeq.ord, env.var = merge_samdat$Precipitation.dormant, plot = FALSE)
  vector_plot1 <- vector_plot1$plot+ xlab("NMDS1")+ylab("NMDS2")+coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))+scale_color_continuous(name="Precipitation (mm)")+xlab("NMDS1")+ylab("NMDS2")+coord_fixed()
  
  bdiv_species<- plot_ordination(ordination = physeq.ord, physeq = top20, type = "species")
  bdiv_species$data$ASV <- rownames(bdiv_species$data)
  species_plot <- bdiv_species+geom_text_repel(aes(label=ASV))
  
  species_key <- species_plot$data
  
  figure<- sample_plot / vector_plot1 # +plot_annotation(tag_levels = "A")
  
  return (
    list(
      figure,
      species_plot,
      species_key,
      results_perm,
      results_envfit
    )
  )
}

physeq_mantel <- function(physeq){
  otutab <- otu_table(physeq) %>% transform(transform="hellinger") %>% as.matrix()
  
  xy <- physeq %>% sample_data%>% as_tibble()%>% select("y"= Latitude, "x"=Longitude) 
  
  xydist <- geodist(xy)
  bcdist <- vegdist(otutab, "euclidean")
  mantel(xydist, bcdist,permutations = 1000, method = "spearman")
}

beta_div_16s <- plot_beta_div(must_16s_merge)

Figure_4_bacterial_beta_diversity <- beta_div_16s[[1]]
Supplemental_figure_3_bacteria_beta_taxa_plot <-  beta_div_16s[[2]]

Figure_4_bacterial_beta_diversity <- Figure_4_bacterial_beta_diversity+plot_annotation(tag_levels = "A")
ggsave("Figures_and_tables_check/Figure_4_bacterial_beta_diversity.pdf", Figure_4_bacterial_beta_diversity, width = 8,height = 8,units = "in" )

physeq_mantel(must_16s_2016)
physeq_mantel(must_16s_2017)

beta_div_its <- plot_beta_div(must_its_merge)

Figure_5_fungal_beta_diversity <- beta_div_its[[1]]
Supplemental_figure_4_fungal_beta_taxa_plot <-  beta_div_its[[2]]

Figure_5_fungal_beta_diversity <- Figure_5_fungal_beta_diversity+plot_annotation(tag_levels = "A")
ggsave("Figures_and_tables_check/Figure_5_fungal_beta_diversity.pdf", Figure_5_fungal_beta_diversity, width = 8,height = 8,units = "in"  )

#Save ASV keys for supplemental Figs 3 and 4
tax_key_16s <- beta_div_16s[[3]]
tax_key_its <- beta_div_its[[3]]
physeq_mantel(must_its_2016)
physeq_mantel(must_its_2017)