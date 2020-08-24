# load packages and data

#Load packages
library(dunn.test)
library(dada2)
library(lawstat)
library(vegan)
library(microbiome)
library(microbiomeSeq)
library(geodist)
library(tidyverse)
library(ggordiplots)
library(corncob)
library(ggmap)
library(ggrepel)
library(shades)
library(patchwork)
library(ggsn)
library(phyloseq)
library(WriteXLS)
library(RColorBrewer)
library(cowplot)
library(knitr)
library(lme4)

# read in list of phyloseq objects
phylist <- readRDS("phyloseq_objects.rds")

# extract phyloseq objects from list
must_16s <- phylist[[1]]
must_16s_merge <- phylist[[2]]
must_its <- phylist[[3]]
must_its_merge <- phylist[[4]]



# convert merged and unmerged sample data to tibbles for easier access. make AVA an ordered factor
samdat <- sample_data(must_16s) %>% as_tibble()
merge_samdat <- read.csv("merge_samdat.csv", row.names = 1)
merge_samdat$AVA <- factor(merge_samdat$AVA, levels=c("Santa Barbara","Monterey", "Sonoma","Mendocino","Willamette Valley"), ordered = TRUE)
merge_samdat$Shipment <- rownames(merge_samdat)
sample_data(must_16s_merge) <- sample_data(merge_samdat)
sample_data(must_its_merge) <- sample_data(merge_samdat)


# some analyses are conducted separately for each year
must_16s_2016 <- subset_samples(must_16s_merge, Year == "2016")
must_16s_2017 <- subset_samples(must_16s_merge, Year == "2017")
must_its_2016 <- subset_samples(must_its_merge, Year == "2016")
must_its_2017 <- subset_samples(must_its_merge, Year == "2017")

### plotting theme
# set theme
theme_set(theme_bw())
# set palette for sites
mypal <- brightness(brewer.pal(5,"Spectral"), 0.8)
# shape values
shapevals <- c(21:25)


gg.gap <- function (plot, ylim, segments, tick_width, rel_heights, vjust = 0, 
          margin = c(top = 0, right = 0, bottom = 0, left = 1), ...) 
{
  if (!is.list(segments)) {
    segments = list(segments)
  }
  if (all(missing(ylim), is.null(plot$coordinates$limits$y))) {
    stop("ylim is undefined")
  }
  else if (ylim[1] == ylim[2]) {
    stop("ylim should not be the same number")
  }
  else if (missing(ylim)) {
    ylim = plot$coordinates$limits$y
  }
  for (j in 1:length(segments)) {
    seg1 = segments[[j]][1]
    seg2 = segments[[j]][2]
    if (seg1 > seg2) {
      if (ylim[1] < ylim[2]) {
        msg = paste0("No.", j, " segment: c(", seg1, 
                     ",", seg2, ") is wrong. It should be ", "c(", 
                     seg2, ",", seg1, ")")
        stop(msg)
      }
    }
    else if (seg1 < seg2) {
      if (ylim[1] > ylim[2]) {
        msg = paste0("No.", j, " segment: c(", seg1, 
                     ",", seg2, ") is wrong. It should be ", "c(", 
                     seg2, ",", seg1, ")")
        stop(msg)
      }
    }
    else if (seg1 == seg2) {
      msg = paste0("No.", j, " segment: c(", seg1, ",", 
                   seg2, ") is wrong. tick_width should not be equal")
      stop(msg)
    }
  }
  if (length(segments) >= 2) {
    if (ylim[1] < ylim[2]) {
      for (k in 2:length(segments)) {
        pre.2 = segments[[k - 1]][2]
        suf.1 = segments[[k]][1]
        if (pre.2 > suf.1) {
          pre = paste0("c(", segments[[k - 1]][1], ",", 
                       segments[[k - 1]][2], ")")
          suf = paste0("c(", segments[[k]][1], ",", segments[[k]][2], 
                       ")")
          msg = paste0("Segments ", k - 1, " and ", k, 
                       ": ", pre, ",", suf, " are wrong. They should be ", 
                       suf, ",", pre)
          stop(msg)
        }
      }
    }
    else if (ylim[1] > ylim[2]) {
      for (k in 2:length(segments)) {
        pre.2 = segments[[k - 1]][2]
        suf.1 = segments[[k]][1]
        if (pre.2 < suf.1) {
          pre = paste0("c(", segments[[k - 1]][1], ",", 
                       segments[[k - 1]][2], ")")
          suf = paste0("c(", segments[[k]][1], ",", segments[[k]][2], 
                       ")")
          msg = paste0("Segments ", k - 1, " and ", k, 
                       ": ", pre, ",", suf, " are wrong. They should be ", 
                       suf, ",", pre)
          stop(msg)
        }
      }
    }
  }
  if (ylim[1] < ylim[2]) {
    if (min(unlist(segments)) <= ylim[1]) 
      stop("the minimum of segments must be more than the minium of ylim")
    if (max(unlist(segments)) > ylim[2]) 
      stop("the maximum of segments must be lower than maximum of ylim")
  }
  else if (ylim[1] > ylim[2]) {
    if (min(unlist(segments)) <= ylim[2]) 
      stop("the minimum of segments must be more than the minium of ylim")
    if (max(unlist(segments)) > ylim[1]) 
      stop("the maximum of segments must be lower than maximum of ylim")
  }
  if (missing(tick_width)) {
    tick_width = rep(abs(ylim[2] - ylim[1])/10, (length(segments) + 
                                                   1))
  }
  if ((length(tick_width) - length(segments)) < 1) {
    int_len = length(tick_width)
    for (m in (int_len + 1):(length(segments) + 1)) {
      tick_width[m] = tick_width[int_len]
    }
  }
  seg_heights = 0
  y_heights = 1
  if (length(seg_heights) < length(segments)) {
    seg_heights_len = length(seg_heights)
    for (m in (seg_heights_len + 1):length(segments)) {
      seg_heights[m] = seg_heights[seg_heights_len]
    }
  }
  if (length(y_heights) < (length(segments) + 1)) {
    y_heights_len = length(y_heights)
    for (m in (y_heights_len + 1):(length(segments) + 1)) {
      y_heights[m] = y_heights[y_heights_len]
    }
  }
  if (length(plot$scales$scales) == 0) {
    trans = "identity"
  }
  else if ("trans" %in% names(plot$scales$scales[[1]])) {
    trans = plot$scales$scales[[1]]$trans
  }
  else {
    trans = "identity"
  }
  if ("reverse" %in% trans) {
    if (ylim[1] < ylim[2]) {
      msg = paste0("ylim: ", "c(", ylim[1], ",", ylim[2], 
                   ")", " is wrong. It should be ", "c(", ylim[2], 
                   ",", ylim[1], ")")
      stop(msg)
    }
  }
  if ("identity" %in% trans) {
    if (ylim[1] > ylim[2]) {
      msg = paste0("ylim: ", "c(", ylim[1], ",", ylim[2], 
                   ")", " is wrong. It should be ", "c(", ylim[2], 
                   ",", ylim[1], ")")
      stop(msg)
    }
  }
  for (i in 1:length(segments)) {
    gap = unlist(segments[i])
    if (i == 1) {
      if (ylim[1] < ylim[2]) {
        breaks = seq(ylim[1], gap[1], by = tick_width[i])
      }
      else if (ylim[1] > ylim[2]) {
        breaks = seq(gap[1], ylim[1], by = tick_width[i])
      }
      p_segment.i <- plot + coord_cartesian(ylim = c(ylim[1], 
                                                     gap[1])) + theme(panel.border = element_blank()) + 
        theme(
          plot.title = element_blank(), legend.position = "none") + 
        scale_y_continuous( trans = trans, 
                            breaks = breaks, labels = scales::number_format(accuracy = 0.01)) + ylab(label = NULL)
      p_segment = list(p_segment.i)
      names(p_segment)[length(p_segment)] = i
      rel_heigh = c(y_heights[i], seg_heights[i])
    }
    else {
      if (ylim[1] < ylim[2]) {
        breaks = seq(ylim[1], gap[1], by = tick_width[i])
      }
      else if (ylim[1] > ylim[2]) {
        breaks = seq(gap[1], ylim[1], by = tick_width[i])
      }
      p_segment.i <- plot + coord_cartesian(ylim = c(unlist(segments[i - 
                                                                       1])[2], gap[1])) + theme(panel.border = element_blank()) + 
        theme(axis.line.y = element_line(), legend.position = "none", 
              axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
              title = element_blank(), axis.title.x = element_blank()) + 
        scale_y_continuous( breaks = breaks, 
                            trans = trans, labels = scales::number_format(accuracy = 0.01)) + ylab(label = NULL)
      p_segment = c(p_segment, list(NULL), list(p_segment.i))
      names(p_segment)[length(p_segment)] = i
      rel_heigh = c(rel_heigh, y_heights[i], seg_heights[i])
    }
    if (i == length(segments)) {
      if (ylim[1] < ylim[2]) {
        breaks = seq(gap[2], ylim[2], by = tick_width[i + 
                                                        1])
      }
      else if (ylim[1] > ylim[2]) {
        breaks = seq(ylim[2], gap[2], by = tick_width[i + 
                                                        1])
      }
      p_segment.i <- plot + coord_cartesian(ylim = c(gap[2], 
                                                     ylim[2])) + theme(panel.border = element_blank()) + 
        theme(
          legend.position = "none", axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
        scale_y_continuous( breaks = breaks, trans = trans, labels = scales::number_format(accuracy = 0.01)) + 
        ylab(label = NULL)
      p_segment = c(p_segment, list(NULL), list(p_segment.i))
      names(p_segment)[length(p_segment)] = i + 1
      rel_heigh = c(rel_heigh, y_heights[i])
    }
  }
  p_segment = rev(p_segment)
  if (missing(rel_heights)) {
    rel_heights = rev(rel_heigh)
  }
  else {
    rel_heights = rev(rel_heights)
  }
  if (is.null(plot$theme$axis.title.y$angle)) {
    angle = 90
  }
  else {
    angle = plot$theme$axis.title.y$angle
  }
  plot_grid(plotlist = p_segment, ncol = 1, align = "v", rel_heights = rel_heights) + 
    theme(plot.margin = unit(margin, "cm")) + draw_label(label = plot$labels$y, 
                                                         x = 0, hjust = plot$theme$axis.title.y$hjust, vjust = vjust, 
                                                         fontfamily = plot$theme$axis.title.y$family, fontface = plot$theme$axis.title.y$face, 
                                                         size = 9, angle = angle, lineheight = plot$theme$axis.title.y$lineheight, 
                                                         colour = plot$theme$axis.title.y$colour)
}
