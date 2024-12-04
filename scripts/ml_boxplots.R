# load libraries
library(tidyverse)
library(ggpubr)
library(ggplot2)

setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/machine_learning/")
disease.all = read.csv("diseased-all_filtmags/diseased-all_performance_results.csv")
disease.all$group = "Diseased (all)"
disease.inf = read.csv("diseased-inf_filtmags/diseased-inf_performance_results.csv")
disease.inf$group = "Diseased (infection only)"
ml.data = rbind(disease.all, disease.inf)
ml.melt = reshape2::melt(ml.data)
ml.melt = ml.melt[which(ml.melt$variable %in% c("AUC", "prAUC")),]
ml.melt$variable = gsub("AUC", "AUROC", ml.melt$variable)
ml.melt$variable = gsub("prAUROC", "AUPRC", ml.melt$variable)
ml.melt$variable = factor(ml.melt$variable, levels=c("AUROC", "AUPRC"))

# rename methods
ml.melt$method = recode(ml.melt$method, "glmnet" = "Ridge Regression", "rf" = "Random Forest", "xgbTree" = "Gradient Boosting")

# AUC scores for different methods
ml.plot = ggplot(ml.melt, aes(x = method, y = value, fill=group)) +
  geom_point(alpha = 0.8, color = "darkgrey", size = 1, position = position_jitterdodge(jitter.width=0.1, dodge.width=0.5, jitter.height = 0)) +
  geom_boxplot(alpha=0.5, outlier.colour=NA, width=0.5) + # This removes the outliers
  geom_hline(yintercept = 0.7, linetype="dashed") +
  coord_flip() +
  facet_wrap(~ variable) +
  scale_fill_manual(values = c("steelblue3", "tomato3")) +
  theme_minimal() +
  ylim(0.5,1) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=14),
    axis.line = element_line(linewidth=0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size=12),
    panel.spacing = unit(2, "lines"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14, margin = margin(r = 20)),
    axis.text = element_text(size = 12),
    axis.ticks = element_line(color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  ylab("Score")
ggsave(file="../figures/ml_box.pdf", height=6, width=9)
