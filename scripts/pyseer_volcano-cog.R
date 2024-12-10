# load libraries
library(ggplot2)
library(reshape2)
library(ggrepel)
library(grid)
library(tidyverse)
library(ggpubr)

setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/gwas/")
COG_pyseer_health = read.csv("diseased-inf_results-cogs.csv")
pyseer_mags = read.delim("gtype_results-genes.tsv")
pyseer_mags$p.adjust = p.adjust(pyseer_mags$lrt.pvalue, method="fdr")
pyseer_mags = pyseer_mags[which(pyseer_mags$p.adjust < 0.05),]
pyseer_health = read.delim("diseased-inf_results-genes.tsv")
pyseer_health$p.adjust = p.adjust(pyseer_health$lrt.pvalue, method="fdr")
pyseer_health = pyseer_health[which(!pyseer_health$variant %in% pyseer_mags$variant),]
pyseer_health$Direction = ifelse(pyseer_health$p.adjust < 0.05, ifelse(pyseer_health$beta > 0, "Diseased", "Healthy"), "Not Sig")
sig.genes = pyseer_health[which(pyseer_health$Direction != "Not Sig"),"variant"]

#function to create volcano plot
volcano = ggplot(pyseer_health, aes(x = beta, y = -log10(p.adjust), color = Direction)) + 
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Diseased" = "#D55E00", "Healthy" = "#0072B2")) +
  ylab(bquote("-log"[10]*"(FDR)")) +
  xlab("Effect size (beta)") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(linewidth=0.5),
    plot.title = element_blank(),
    legend.title = element_blank(),
    
    # Show legend title
    legend.text = element_text(size = 16),
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "white")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

# cog analysis
COG_descriptions = c(
  B = "Chromatin Structure and dynamics",
  J = "Translation",
  L = "Replication, recombination & repair",
  K = "Transcription",
  O = "Molecular chaperones and related functions",
  M = "Cell wall structure & outer membrane",
  N = "Secretion, motility and chemotaxis",
  T = "Signal transduction",
  P = "Inorganic ion transport and metabolism",
  U = "Intracellular trafficking and secretion",
  C = "Energy production and conversion",
  G = "Carbohydrate metabolism and transport",
  E = "Amino acid metabolism and transport",
  F = "Nucleotide metabolism and transport",
  H = "Coenzyme metabolism",
  I = "Lipid metabolism",
  D = "Cell division and chromosome partitioning",
  R = "General functional prediction only",
  S = "No functional prediction",
  Q = "Secondary Structure",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  Z = "Cytoskeleton"
)

# transform COGs using descriptions
transform_cogs = function(df, descriptions) {
  df$COG = as.character(df$COG)
  df$COG = sapply(df$COG, function(x) if (x %in% names(descriptions)) return(descriptions[x]) else return(x))
  return(df)
}

COG_pyseer_health = transform_cogs(COG_pyseer_health, COG_descriptions)
COG_pyseer_health$Positive = COG_pyseer_health$Positive/sum(COG_pyseer_health$Positive)*100
COG_pyseer_health$Negative = COG_pyseer_health$Negative/sum(COG_pyseer_health$Negative)*100
COG_pyseer_health = COG_pyseer_health[which(COG_pyseer_health$COG != "No functional prediction"),]

# prepare dataframe for plot
data_long = melt(COG_pyseer_health, id.vars = "COG", variable.name = "Type", value.name = "Count")
cog_sum = aggregate(Count ~ COG + Type, data = data_long, sum)
cog_sum = reshape2::dcast(cog_sum, COG ~ Type, value.var = "Count")
cog_sum$Sum = cog_sum$Positive + cog_sum$Negative
cog_order = cog_sum[order(cog_sum$Sum), "COG"]
data_long$COG = factor(data_long$COG, levels = cog_order)

# plot data
cog.plot = ggplot(data_long, aes(x = COG, y = Count, fill = Type)) +
  geom_bar(data = subset(data_long, Type == "Positive"), stat = "identity") +
  geom_bar(data = subset(data_long, Type == "Negative"), stat = "identity", aes(y = -Count)) +
  scale_y_continuous(labels = abs) +
  coord_flip() +
  ylab("Proportion of significant genes (%)") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=16),
    axis.line = element_line(linewidth=0.5),
    legend.text = element_text(size=16),
    legend.title = element_blank(),
    legend.position = "bottom") +
  scale_fill_manual(values = c("Positive" = "#D55E00", "Negative" = "#0072B2"), labels = c("Diseased", "Healthy"), name = "")

# combine
gwas.fig = ggarrange(volcano, cog.plot, common.legend=TRUE, legend = "bottom", widths=c(1,1.2), labels=c("a", "b"), font.label = list(size=18))