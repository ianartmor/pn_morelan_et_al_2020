library(corncob)
library(phyloseq)
library(purrr)
library(dplyr)
library(forcats)
data("soil_phylo")


soil_diftest_res <- differentialTest(
  data = soil_phylo %>% tax_glom("Genus"), formula = ~DayAmdmt, phi.formula = ~1,
  formula_null = ~1, phi.formula_null = ~1, test = "LRT", full_output = TRUE
)

differentialTestPredictions <- function(diftest_res, B = 1000, taxa) {
  
  # Below line is probably not ideal! Would be safer if differentialTest outputs were named lists
  names(diftest_res$full_output) <- taxa_names(diftest_res$data)
  
  # default is all taxa
  if (missing(taxa)) {
    taxa <- taxa_names(diftest_res$data)
  } else {
    taxa <- taxa
  }
  
  purrr::map(
    diftest_res$full_output[taxa],
    ~ corncob:::plot.bbdml(.,
                           B = B, data_only = TRUE)
  )
}

#low bootstrap
prediction_list <- differentialTestPredictions(
  diftest_res = soil_diftest_res,
  B = 5,
  taxa = soil_diftest_res$significant_taxa
)

## plot
samdat <- sample_data(soil_phylo) %>%
  data.frame() %>%
  rownames_to_column("samples")

prediction_list_with_samdat <- prediction_list %>% map(~ left_join(., samdat, by = "samples"))

plotlist <- prediction_list_with_samdat %>%
  map2(.x = ., .y=names(.),
    ~ mutate(.data = ., samples = fct_reorder(samples, DayAmdmt)) %>%
      ggplot(data = .) +
      geom_pointrange(aes(x = samples, ymin = ymin, ymax = ymax, y = RA, color = DayAmdmt)) +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ggtitle(.y)
  )

plotlist[1]
