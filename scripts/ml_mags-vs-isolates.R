# load libraries
library(tidyverse)
library(data.table)
library(ggsignif)
library(ggplot2)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/machine_learning")

input.files = list.files(path = ".", pattern = "_results.csv", recursive = FALSE)
load_ml = function(x) {
  base = basename(x)
  type = gsub("_results.csv", "", base)
  if(grepl("isolates", type)) {
    label = "Isolates"
    type = gsub("_isolates", "", type)
  } else {
    label = "MAGs + Isolates"
  }
  ml.in = read.csv(x)
  ml.in$Data = label
  ml.in$Disease = type
  return(ml.in)
}

ml.list = lapply(input.files, load_ml)
ml.combined = as.data.frame(rbindlist(ml.list))
ml.melt = reshape2::melt(ml.combined)

# plot
ml.subset = ml.melt[which(ml.melt$variable %in% c("AUC", "prAUC") & ml.melt$Disease == "diseased-inf"),]
ml.subset$variable = recode(ml.subset$variable, "AUC" = "AUROC", "prAUC" = "AUPRC")
box.plot = ggplot(ml.subset, aes(x=Data, y=value, fill=Data)) +
  geom_boxplot(alpha=0.5, outlier.shape=NA, width=0.6) +
  geom_point(alpha=0.8, size=0.5, colour="grey40", position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0)) +
  geom_hline(yintercept=0.7, linetype = "dashed") +
  scale_fill_manual(values=rev(c("darkolivegreen3", "tomato3"))) +
  facet_wrap(~ variable) +
  coord_cartesian(ylim=c(0.4, 1.05), expand = TRUE) +
  ylab("Score") +
  theme_classic() +
  theme(panel.spacing = unit(2, "lines")) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(size=16)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=14)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(
    comparisons = list(c("MAGs + Isolates", "Isolates")),
    map_signif_level = FALSE, textsize=4.5
  )

# combine with roc curve
source("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/scripts/alex/ml_roc-curves.R")
ml.fig = ggarrange(roc.curve, box.plot, widths=c(1,1.1), labels=c("c", "d"), font.label=list(size=18))
