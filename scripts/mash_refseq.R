# load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/mash/")
metadata = read.delim("../metadata/Metadata_15022025.tsv")
mags = metadata$Genome[which(metadata$Genome_Type == "MAG")]
rownames(metadata) = metadata$Genome
metadata[is.na(metadata)] = "Unknown"
distances = fread("distances.tab")
str(distances)
setnames(distances, c("ref_genome", "mag", "distance", "p_value", "shared_hashes"))

# get top hits only
closest_matches = distances %>%
  group_by(mag) %>%
  slice_min(order_by = distance, with_ties = FALSE) %>%
  ungroup()

# simplify names
extract_up_to_second_underscore = function(string) {
  parts = strsplit(string, "_")[[1]]
  paste(parts[1:2], collapse = "_")
}
closest_matches$ref_genome = as.vector(sapply(basename(closest_matches$ref_genome), extract_up_to_second_underscore))
closest_matches$mag = as.vector(gsub(".fa.gz", "", basename(closest_matches$mag)))
closest_matches = closest_matches[which(closest_matches$mag %in% mags),]

mash.plot = ggplot(closest_matches, aes(x = distance)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "grey", color = "black", alpha = 0.7) +
  geom_density(color = "#377EB8", linewidth = 1) +
  geom_vline(xintercept = 0.005, linetype="dashed") +
  labs(x = "Mash distance", y = "Number of MAGs") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)))

# assess new lineages
new.lineages = closest_matches[which(closest_matches$distance > 0.005),"mag"]$mag
new.lineages.meta = metadata[new.lineages,c("Country", "Health_Status")]
new.lineages.melt = reshape2::melt(as.matrix(new.lineages.meta))
new.lineages.melt$Count = 1
new.lineages.melt = aggregate(Count ~ value + Var2, data=new.lineages.melt, FUN=sum)
new.lineages.melt$value = factor(new.lineages.melt$value, levels=new.lineages.melt[order(new.lineages.melt$Count),"value"])

# color dictionary
countries = rev(brewer.pal(10, "Paired"))
names(countries) = as.vector(new.lineages.melt[which(new.lineages.melt$Var2 == "Country"),"value"])
health = c("#cd5f08", "#0072B2", "grey")
names(health) = c("Diseased", "Healthy", "Unknown")
colors.dict = c(countries, health)

# plot metadata distribution
meta.plot = ggplot(new.lineages.melt, aes(x=reorder(Var2, -Count), y=Count, fill=value)) +
  geom_bar(stat="identity", alpha=0.7) +
  theme_classic() +
  scale_fill_manual(values=colors.dict, name="Variable") +
  scale_x_discrete(labels=c("Health status", "Country")) +
  ylab("Number of MAGs") +
  xlab("") +
  theme(
    axis.title = element_text(size=14),
    axis.text = element_text(size=12)
  )

# combine plot and save
comb.plot = ggarrange(mash.plot, meta.plot, ncol=2, widths=c(1.7,1), labels=c("a", "b"), font.label = list(size=16))
ggsave("../figures/mash_comb.pdf", height=5, width=10)

# save table of best hits
write.table(closest_matches, file="best_matches.tsv", row.names=FALSE, quote=FALSE, sep="\t")