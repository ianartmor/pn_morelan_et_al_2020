# library(dada2)
# library(microbiomeSeq)
# library(microbiome)
# library(corncob)
# library(patchwork)
# library(cowplot)
# library(tidyverse)
# library(phyloseq)
#source("Analysis_scripts/00_setup.R")
#source("Analysis_scripts/07a_da_functions.R")


# differential abundance using corncob

set.seed(5)


theme_set(theme_bw())
merge_gen_16s <- taxa_level(must_16s_merge, "Genus")
toptax <- prune_taxa(top_taxa(merge_gen_16s, 12), merge_gen_16s)

# outlier_detect(toptax)

# no_outliers_16s <- toptax %>%
#  subset_samples(Shipment!="RRV1_2016") %>%
#  subset_samples(Shipment!="SMV1_2017") %>%
#  subset_samples(Shipment != "SMV2_2017") %>%
#  subset_samples(Shipment != "RRV3_2016")

### VINTAGE

vintage_res_16s <- vintage_dif(toptax)

sig_vintage_16s <- toptax %>%
  transform("compositional") %>%
  prune_taxa(taxa = vintage_res_16s$significant_taxa, x = .)


# extract coefficients
vintage_coefficients_16s <- vintage_res_16s$significant_models %>%
  purrr::map(coef) %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  purrr::map_df(., .f = ~ filter(., rowname == "mu.Year")) %>%
  mutate(rowname = vintage_res_16s$significant_taxa)

# toptax and plot
vintage_plot_16s <- sig_vintage_16s %>%
  psmelt() %>%
  filter(!grepl("f__", OTU)) %>%
  split(.$OTU) %>%
  purrr::map2(.x = ., .y = names(.), .f = ~ ggplot(., aes(x = as.factor(Year), y = Abundance * 100, group = Year, fill = AVA, shape = as.factor(Year))) +
    geom_point(show.legend = FALSE) +
    geom_boxplot(alpha = 0.2, show.legend = FALSE) +
    xlab("") +
    ylab("") +
    ggtitle({
      .y
    }) +
    theme(
      plot.title = element_text(size = 8, face = "italic"),
      axis.title = element_text(size = 8, face = "plain"),
      panel.border = element_blank(),
      # legend.title = element_text(size=9),
      # legend.text = element_text(size=7), legend.spacing.y = unit(1, "mm"), legend.box = "horizontal"
    ) +
    scale_fill_manual(values = mypal) +
    scale_shape_manual(values = shapevals, name = "Vintage") +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_fill_manual(values = mypal))

# add URA label to only leftmost plot

vintage_plot_16s[[1]] <- vintage_plot_16s[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))


# use modified gg.gap() function (see setup.R for code) to make gapped plots so outliers can be shown on same graph.
# rel heights parameter used to try to get ticks to be approximately the same spacing. 
vintage_plot_16s[[1]] <- gg.gap( vintage_plot_16s[[1]], ylim = c(0,16), segments = c(4,12), tick_width = 2, rel_heights = c(1,0,0.9))
vintage_plot_16s[[1]]<- ggdraw(vintage_plot_16s[[1]])+draw_label("//", x=0.52,y=0.54, angle=305)

vintage_plot_16s[[2]] <- gg.gap( vintage_plot_16s[[2]], ylim = c(0,19), segments = c(5,18), rel_heights = c(2.8,0,1), tick_width = 1,margin = c(top = 0, right = 0, bottom = 0, left = 0))
vintage_plot_16s[[2]]<- ggdraw(vintage_plot_16s[[2]])+draw_label("//", x=0.33,y=0.74, angle=305)

# rel width needs to be adjusted slightly due to leftmost plot having a y axis label
vintage_plots_16s <- plot_grid(plotlist = vintage_plot_16s, nrow=1, rel_widths = c(1.2,1.1,1,1,1,1))


#vintage_plots_16s <- plot_grid(plotlist = vintage_plot_16s, nrow = 1)

### PRECIPITATION
precip_res_16s <- must_precip_growing_dif(toptax)

# toptax
must_precip_sig_16s <- toptax %>%
  transform("compositional") %>%
  prune_taxa(precip_res_16s$significant_taxa, .)

# split and plot
precip_plot_16s <- must_precip_sig_16s %>%
  psmelt() %>%
  filter(!grepl("f__", OTU)) %>%
  split(.$OTU) %>%
  purrr::map2(.x = ., .y = names(.), .f = ~ ggplot(., aes(x = Precipitation.growing, y = Abundance * 100, fill = AVA, shape = as.factor(Year))) +
    geom_point() +
    # geom_smooth(method="lm", se=FALSE)+
    ggtitle({
      .y
    }) +
    xlab("Precipitation (mm)") +
    ylab("") +
    theme(
      plot.title = element_text(size = 9, face = "italic"),
      panel.border = element_blank(),
      axis.title = element_text(size = 9, face = "plain"), legend.position = "none"
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_fill_manual(values = mypal) +
    scale_shape_manual(values = shapevals))
precip_plots_16s <- plot_grid(plotlist = precip_plot_16s, nrow = 1)

precip_coefficients_16s <- precip_res_16s$significant_models %>%
  purrr::map(coef) %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  purrr::map_df(., .f = ~ filter(., rowname == "mu.Precipitation.growing")) %>%
  mutate(rowname = precip_res_16s$significant_taxa)


### TA
ta_res_16s <- must_ta_dif(toptax)

# toptax
must_ta_sig_16s <- toptax %>%
  transform("compositional") %>%
  prune_taxa(ta_res_16s$significant_taxa, .)

# split and plot
ta_plot_16s <- must_ta_sig_16s %>%
  psmelt() %>%
  filter(!grepl("f__", OTU)) %>%
  split(.$OTU) %>%
  purrr::map2(.x = ., .y = names(.), .f = ~ ggplot(., aes(x = Must_TA, y = Abundance * 100, fill = AVA, shape = as.factor(Year))) +
    geom_point() +
    # geom_smooth(method="lm", se=FALSE)+
    ggtitle({
      .y
    }) +
    xlab("TA (g/L)") +
    ylab("") +
    theme(
      plot.title = element_text(size = 9, face = "italic"),
      panel.border = element_blank(),
      axis.title = element_text(size = 9, face = "plain"), legend.position = "none"
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_fill_manual(values = mypal) +
    scale_shape_manual(values = shapevals))

ta_plot_16s[[1]] <- ta_plot_16s[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))

ta_plots_16s <- plot_grid(plotlist = ta_plot_16s, nrow = 1)


ta_coefficients_16s <- ta_res_16s$significant_models %>%
  purrr::map(coef) %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  purrr::map_df(., .f = ~ filter(., rowname == "mu.Must_TA")) %>%
  mutate(rowname = ta_res_16s$significant_taxa)

#### PH

ph_res_16s <- must_pH_dif(toptax)

# toptax
must_ph_sig_16s <- toptax %>%
  transform("compositional") %>%
  prune_taxa(ph_res_16s$significant_taxa, .)


# split and plot
ph_plot_16s <- must_ph_sig_16s %>%
  psmelt() %>%
  filter(!grepl("f__", OTU)) %>%
  split(.$OTU) %>%
  map2(.x = ., .y = names(.), .f = ~ ggplot(., aes(x = Must_pH, y = Abundance * 100, fill = AVA, shape = as.factor(Year))) +
    geom_point() +
    # geom_smooth(method="lm", se=FALSE)+
    xlab("pH") +
    ylab("") +
    ggtitle({
      .y
    }) +
    theme(
      plot.title = element_text(size = 9, face = "italic"), axis.title = element_text(size = 9, face = "plain"),
      panel.border = element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    scale_fill_manual(values = mypal) +
    scale_shape_manual(values = shapevals))

ph_plot_16s[[1]] <- ph_plot_16s[[1]]+ ylab("URA (%)")+theme(axis.title = element_text(size=9))


# add gap for outlier
ph_plot_16s[[1]]<- gg.gap(ph_plot_16s[[1]], ylim=c(0,22.25), segments = c(2,21.75), rel_heights = c(2.4,0,1), tick_width = 0.5)
ph_plot_16s[[1]]<- ggdraw(ph_plot_16s[[1]])+draw_label("//", x=0.3,y=0.71, angle=305)




ph_plots_16s <- plot_grid(plotlist = ph_plot_16s, nrow = 1)

# extract coefficients
ph_coefficients_16s <- ph_res_16s$significant_models %>%
  purrr::map(coef) %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  purrr::map_df(., .f = ~ filter(., rowname == "mu.Must_pH")) %>%
  mutate(rowname = ph_res_16s$significant_taxa)


### TSS

tss_res_16s <- must_tss_dif(toptax)
# toptax
must_tss_sig_16s <- toptax %>%
  transform("compositional") %>%
  prune_taxa(tss_res_16s$significant_taxa, .)

# split and plot
tss_plot_16s <- must_tss_sig_16s %>%
  psmelt() %>%
  mutate(Vintage = as.factor(Year)) %>%
  filter(!grepl("f__", OTU)) %>%
  split(.$OTU) %>%
  map2(.x = ., .y = names(.), .f = ~ ggplot(., aes(x = Must_TSS, y = Abundance * 100, fill = AVA, shape = Vintage)) +
    geom_point() +
    # geom_smooth(method="lm", se=FALSE)+
    xlab("TSS (°Brix)") +
    ylab("") +
    ggtitle({
      .y
    }) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
    theme(
      plot.title = element_text(size = 9, face = "italic"),
      panel.border = element_blank(), legend.position = "none",
      axis.title = element_text(size = 9, face = "plain")
    ) +
    scale_fill_manual(values = mypal) +
    scale_shape_manual(values = shapevals))

# add y axis label
tss_plot_16s[[1]] <- tss_plot_16s[[1]]+ylab("URA (%)")+theme(axis.title = element_text(size=9))

# add gap for outlier
tss_plot_16s[[1]]<- gg.gap(tss_plot_16s[[1]], ylim=c(0,22.25), segments = c(2,21.75), rel_heights = c(2.4,0,1), tick_width = 0.5)
tss_plot_16s[[1]]<- ggdraw(tss_plot_16s[[1]])+draw_label("//", x=0.20,y=0.71, angle=305)
tss_plots_16s <- plot_grid(plotlist = tss_plot_16s, nrow = 1)

# extract coefficients
tss_coefficients_16s <- tss_res_16s$significant_models %>%
  purrr::map(coef) %>%
  map(data.frame) %>%
  map(rownames_to_column) %>%
  purrr::map_df(., .f = ~ filter(., rowname == "mu.Must_TSS")) %>%
  mutate(rowname = tss_res_16s$significant_taxa)

# combine plots
bact_difabund <- (vintage_plots_16s / precip_plots_16s / ta_plots_16s / ph_plots_16s / tss_plots_16s) + plot_annotation(tag_levels = "A")



WriteXLS(x = c("vintage_coefficients_16s", "precip_coefficients_16s", "ta_coefficients_16s", "ph_coefficients_16s", "tss_coefficients_16s"), ExcelFileName = "Table_S6_corncob_bact.xls")
ggsave("Figures_and_tables_check/Figure_6_bacteria_differential_abundance.pdf", bact_difabund, width = 8.5, height = 11, units = "in")

