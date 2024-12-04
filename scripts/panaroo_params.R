# load libraries
library(ggplot2)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/")
pan.genomes = read.delim("panaroo/pan-size_params.tsv", header=FALSE)
colnames(pan.genomes) = c("Mode", "Identity", "Paralog", "Pan-genome", "Core-genome")
pan.melt = reshape2::melt(pan.genomes, id.vars=c("Mode", "Identity", "Paralog"))
pan.melt$Mode = gsub("mod", "Moderate", pan.melt$Mode)
pan.melt$Mode = gsub("str", "Strict", pan.melt$Mode)
pan.melt$Paralog = gsub("NMP", "Not-merged", pan.melt$Paralog)
pan.melt$Paralog = gsub("MP", "Merged", pan.melt$Paralog)

# calculate max variation
max_variation_percentage = function(x) {
  mean_value = mean(x)
  max_deviation = max(abs(x - mean_value))
  (max_deviation / mean_value) * 100
}

pan.max = max_variation_percentage(pan.genomes$`Pan-genome`)
core.max = max_variation_percentage(pan.genomes$`Core-genome`)

# scatterplot
scatter = ggplot(pan.melt, aes(x=Identity, y=value, colour=Mode, shape=Paralog)) +
  geom_point(size=5, alpha=0.5) +
  geom_line() +
  facet_wrap(~ variable, scales="free") +
  scale_x_continuous(breaks=c(90,95), expand=c(0.5,0.5)) +
  ylab("Number of genes") +
  xlab("Identity threshold (%)") +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth=0.5),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14),
    strip.text = element_text(size=14),
    legend.text = element_text(size=10),
    legend.title = element_text(size=10),
    panel.spacing = unit(1, "cm"))
ggsave(file = "figures/panaroo_params.pdf", height=5, width=9)


                                 
