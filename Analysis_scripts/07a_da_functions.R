
### DIFFERENTIAL ABUNDANCE MULTIPLE TAXA
# differential abundance across vintages

vintage_dif <- function(ps) {
  differentialTest(
    formula = ~ Year + AVA,
    phi.formula = ~AVA,
    formula_null = ~AVA,
    phi.formula_null = ~AVA,
    data = ps,
    test = "LRT", boot = FALSE,
    full_output = TRUE,
    fdr_cutoff = 0.05
  )
}

# differential abundance with precipitation, controlling for AVA and vintage
must_precip_growing_dif <- function(ps) {
  differentialTest(
    formula = ~ Precipitation.growing + AVA + Year,
    phi.formula = ~ AVA + Year,
    formula_null = ~ AVA + Year,
    phi.formula_null = ~ AVA + Year,
    data = ps,
    test = "LRT", boot = FALSE,
    full_output = TRUE,
    fdr_cutoff = 0.05
  )
}

# differential abundance with titratable acidity, controlling for AVA and vintage
must_ta_dif <- function(ps) {
  differentialTest(
    formula = ~ Must_TA + Year + AVA,
    phi.formula = ~ AVA + Year,
    formula_null = ~ AVA + Year,
    phi.formula_null = ~ AVA + Year,
    data = ps,
    test = "LRT", boot = FALSE,
    full_output = TRUE,
    fdr_cutoff = 0.05
  )
}

# differential abundance with pH, controlling for AVA and vintage
must_pH_dif <- function(ps) {
  differentialTest(
    formula = ~ Must_pH + Year + AVA,
    phi.formula = ~ AVA + Year,
    formula_null = ~ AVA + Year,
    phi.formula_null = ~ AVA + Year,
    data = ps,
    test = "LRT", boot = FALSE,
    full_output = TRUE,
    fdr_cutoff = 0.05
  )
}


# differential abundance with total soluble solids, controlling for AVA and vintage
must_tss_dif <- function(ps) {
  differentialTest(
    formula = ~ Must_TSS + Year + AVA,
    phi.formula = ~ AVA + Year,
    formula_null = ~ AVA + Year,
    phi.formula_null = ~ AVA + Year,
    data = ps,
    test = "LRT", boot = FALSE,
    full_output = TRUE,
    fdr_cutoff = 0.05
  )
}


### SINGLE TAXON MODELS 

## Vintage
# plots three single taxon models:
# 1: taxon only
# 2: AVA
# 3: AVA + Year with AVA as null
vintage_dif_single <- function(ps, ASV) {
  theme_set(theme_bw())
  asv <- ASV
  taxonly <- deparse(substitute(1))
  ava <- deparse(substitute(AVA))
  full <- deparse(substitute(AVA + Year))

  tax_only_form <- reformulate(taxonly, response = asv)
  ava_form <- reformulate(ava, response = asv)
  full_form <- reformulate(full, response = asv)

  tax_only <- bbdml(
    formula = tax_only_form,
    phi.formula = ~1,
    data = ps
  )

  tax_null <- bbdml(
    formula = ava_form,
    phi.formula = ~AVA,
    data = ps
  )

  tax_full <- bbdml(
    formula = full_form,
    phi.formula = ~ Year + AVA,
    data = ps
  )

  pval <- signif(lrtest(mod_null = tax_null, mod = tax_full), 3)


  p1 <- plot(tax_only,
    color = "Year",
    sample_names = FALSE, B = 1000, data_only = TRUE
  ) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Year)) %>%
    mutate(Year = as.factor(Year)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Year), size = 0.6, fatten = 0.3) +
    ggtitle(paste(ASV, " only", sep = "")) +
    theme(
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_discrete(name = "Vintage")

  p2 <- plot(tax_null, color = "Vintage", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Year)) %>%
    mutate(Year = as.factor(Year)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Year), size = 0.6, fatten = 0.3) +
    ggtitle("AVA as covariate \n(null model)") +
    theme(
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("Shipment (ordered by Vintage)") +
    scale_color_discrete(name = "Vintage")

  p3 <- plot(tax_full, color = "Vintage", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Year)) %>%
    mutate(Year = as.factor(Year)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Year), size = 0.6, fatten = 0.3) +
    ggtitle(paste("AVA+Vintage as covariates, \n p=", as.character(pval), sep = "")) +
    theme(
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_discrete(name = "Vintage")


  multi <- p1 + p2 + p3 + plot_layout(guides = "collect")
  return(multi)
}


## Precipitation
# Plots 3 models:
# 1: Taxon only
# 2: AVA+Year
# 3: AVA+Year+Precipitation with AVA+Year as null

precip_dif_single <- function(ps, ASV) {
  theme_set(theme_bw())
  asv <- ASV
  taxonly <- deparse(substitute(1))
  null <- deparse(substitute(AVA + Year))
  full <- deparse(substitute(AVA + Year + Precipitation.growing))

  tax_only_form <- reformulate(taxonly, response = asv)
  null_form <- reformulate(null, response = asv)
  full_form <- reformulate(full, response = asv)

  tax_only <- bbdml(
    formula = tax_only_form,
    phi.formula = ~1,
    data = ps
  )

  tax_null <- bbdml(
    formula = null_form,
    phi.formula = ~ AVA + Year,
    data = ps
  )

  tax_full <- bbdml(
    formula = full_form,
    phi.formula = ~ Year + AVA + Precipitation.growing,
    data = ps
  )

  pval <- signif(lrtest(mod_null = tax_null, mod = tax_full), 3)

  p1 <- plot(tax_only,
    color = "Precipitation.growing",
    sample_names = FALSE, B = 1000, data_only = TRUE
  ) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Precipitation.growing)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Precipitation.growing), size = 0.6, fatten = 0.3) +
    ggtitle(paste(ASV, " only", sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = "Precipitation (mm)")

  p2 <- plot(tax_null, color = "Precipitation.growing", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Precipitation.growing)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Precipitation.growing), size = 0.6, fatten = 0.3) +
    ggtitle("AVA+Vintage as covariates \n(null model)") +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("Shipment (ordered by precipitation)") +
    scale_color_continuous(name = "Precipitation (mm)")

  p3 <- plot(tax_full, color = "Precipitation.growing", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Precipitation.growing)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Precipitation.growing), size = 0.6, fatten = 0.3) +
    ggtitle(paste("AVA+Vintage+Precipitation as covariates, \n p=", as.character(pval), sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = "Precipitation (mm)")

  multi <- p1 + p2 + p3 + plot_layout(guides = "collect")
  return(multi)
}

# Plots 3 single taxon models:
# 1: Taxon only
# 2: AVA+Year
# 3: AVA+Year+tA with AVA+Year as null

ta_dif_single <- function(ps, ASV) {
  theme_set(theme_bw())
  asv <- ASV
  taxonly <- deparse(substitute(1))
  null <- deparse(substitute(AVA + Year))
  full <- deparse(substitute(AVA + Year + Must_TA))

  tax_only_form <- reformulate(taxonly, response = asv)
  null_form <- reformulate(null, response = asv)
  full_form <- reformulate(full, response = asv)

  tax_only <- bbdml(
    formula = tax_only_form,
    phi.formula = ~1,
    data = ps
  )

  tax_null <- bbdml(
    formula = null_form,
    phi.formula = ~ AVA + Year,
    data = ps
  )

  tax_full <- bbdml(
    formula = full_form,
    phi.formula = ~ Year + AVA + Must_TA,
    data = ps
  )

  pval <- signif(lrtest(mod_null = tax_null, mod = tax_full), 3)

  p1 <- plot(tax_only,
    color = "Must_TA",
    sample_names = FALSE, B = 1000, data_only = TRUE
  ) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_TA)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_TA), size = 0.6, fatten = 0.3) +
    ggtitle(paste(ASV, " only", sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = expression(TA ~ (g ~ tartaric ~ acid ~ L^{
      "-1"
    })))

  p2 <- plot(tax_null, color = "Must_TA", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_TA)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_TA), size = 0.6, fatten = 0.3) +
    ggtitle("AVA+Vintage as covariates \n(null model)") +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("Shipment (ordered by TA conc.)") +
    scale_color_continuous(name = expression(TA ~ (g ~ tartaric ~ acid ~ L^{
      "-1"
    })))

  p3 <- plot(tax_full, color = "Must_TA", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_TA)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_TA), size = 0.6, fatten = 0.3) +
    ggtitle(paste("AVA+Vintage+TA as covariates, \n p=", as.character(pval), sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = expression(TA ~ (g ~ tartaric ~ acid ~ L^{
      "-1"
    })))

  multi <- p1 + p2 + p3 + plot_layout(guides = "collect")
  return(multi)
}

# Plots 3 single-taxon models:
# 1: Taxon only
# 2: AVA+Year
# 3: AVA+Year+pH with AVA+Year as null

pH_dif_single <- function(ps, ASV) {
  theme_set(theme_bw())
  asv <- ASV
  taxonly <- deparse(substitute(1))
  null <- deparse(substitute(AVA + Year))
  full <- deparse(substitute(AVA + Year + Must_pH))

  tax_only_form <- reformulate(taxonly, response = asv)
  null_form <- reformulate(null, response = asv)
  full_form <- reformulate(full, response = asv)

  tax_only <- bbdml(
    formula = tax_only_form,
    phi.formula = ~1,
    data = ps
  )

  tax_null <- bbdml(
    formula = null_form,
    phi.formula = ~ AVA + Year,
    data = ps
  )

  tax_full <- bbdml(
    formula = full_form,
    phi.formula = ~ Year + AVA + Must_pH,
    data = ps
  )

  pval <- signif(lrtest(mod_null = tax_null, mod = tax_full), 3)

  p1 <- plot(tax_only,
    color = "Must_pH",
    sample_names = FALSE, B = 1000, data_only = TRUE
  ) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_pH)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_pH), size = 0.6, fatten = 0.3) +
    ggtitle(paste(ASV, " only", sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = "pH")

  p2 <- plot(tax_null, color = "Must_pH", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_pH)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_pH), size = 0.6, fatten = 0.3) +
    ggtitle("AVA+Vintage as covariates \n(null model)") +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("Shipment (ordered by pH)") +
    scale_color_continuous(name = "pH")

  p3 <- plot(tax_full, color = "Must_pH", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_pH)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_pH), size = 0.6, fatten = 0.3) +
    ggtitle(paste("AVA+Vintage+pH as covariates, \n p=", as.character(pval), sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = "pH")

  multi <- p1 + p2 + p3 + plot_layout(guides = "collect")
  return(multi)
}

# Plots 3 single-taxon models:
# 1: Taxon only
# 2: AVA+Year
# 3: AVA+Year+TSS with AVA+Year as null

tss_dif_single <- function(ps, ASV) {
  theme_set(theme_bw())
  asv <- ASV
  taxonly <- deparse(substitute(1))
  null <- deparse(substitute(AVA + Year))
  full <- deparse(substitute(AVA + Year + Must_TSS))

  tax_only_form <- reformulate(taxonly, response = asv)
  null_form <- reformulate(null, response = asv)
  full_form <- reformulate(full, response = asv)

  tax_only <- bbdml(
    formula = tax_only_form,
    phi.formula = ~1,
    data = ps
  )

  tax_null <- bbdml(
    formula = null_form,
    phi.formula = ~ AVA + Year,
    data = ps
  )

  tax_full <- bbdml(
    formula = full_form,
    phi.formula = ~ Year + AVA + Must_TSS,
    data = ps
  )

  pval <- signif(lrtest(mod_null = tax_null, mod = tax_full), 3)

  p1 <- plot(tax_only,
    color = "Must_TSS",
    sample_names = FALSE, B = 1000, data_only = TRUE
  ) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_TSS)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_TSS), size = 0.6, fatten = 0.3) +
    ggtitle(paste(ASV, " only", sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = "TSS (°Brix)")

  p2 <- plot(tax_null, color = "Must_TSS", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_TSS)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_TSS), size = 0.6, fatten = 0.3) +
    ggtitle("AVA+Vintage as covariates \n(null model)") +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("Shipment (ordered by TSS)") +
    scale_color_continuous(name = "TSS (°Brix)")


  p3 <- plot(tax_full, color = "Must_TSS", sample_names = FALSE, B = 1000, data_only = TRUE) %>%
    rename(Shipment = samples) %>%
    left_join(x = ., y = merge_samdat, by = c("Shipment")) %>%
    mutate(Shipment = fct_reorder(Shipment, .x = Must_TSS)) %>%
    ggplot() +
    geom_pointrange(aes(x = Shipment, ymin = ymin, ymax = ymax, y = RA, color = Must_TSS), size = 0.6, fatten = 0.3) +
    ggtitle(paste("AVA+Vintage+TSS as covariates, \n p=", as.character(pval), sep = "")) +
    theme(
      legend.title.align = 0.5,
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      text = element_text(size = 6)
    ) +
    xlab("") +
    scale_color_continuous(name = "TSS (°Brix)")

  multi <- p1 + p2 + p3 + plot_layout(guides = "collect")
  return(multi)
}


# OUTLIER DETECTION
# returns a data frame with similar dimensions to OTU table with TRUE/FALSE values;
# TRUE where "high" outliers (defined as 3 SDs from mean) were detected. Easily modified to detect low outliers.

outlier_detect <- function(ps) {
  otutab <- ps %>%
    otu_table() %>%
    data.frame()
  
  detect_out <- function(x) {
    highouts <- ((mean(x) + 3 * sd(x)) - x) < 1
    lowouts <- ((mean(x) - 3 * sd(x)) - x) > 1
    return(highouts)
  }
  
  out_tab <- otutab %>%
    mutate(across(.f = detect_out))
  rownames(out_tab) <- rownames(otutab)
  
  return(out_tab)
}

