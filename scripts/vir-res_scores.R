# load libraries
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(ggsignif)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/")
metadata = read.delim("metadata/Metadata_15022025.tsv")
metadata = metadata[which(!is.na(metadata$Health_Status)),]
klebo = read.delim("kleborate/Kleborate_results.tsv")

# only healthy and infection
metadata = metadata[which(metadata$Disease_Name %in% c("Infection", "Diarrhoea", "Healthy")),]

virulence_plot <- ggplot(metadata, aes(x = factor(Virulence_score), fill = Health_Status)) +
  geom_histogram(stat="count", colour="grey", width=0.75, alpha=0.9) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +  
  labs(x = "Virulence score", y = "Number of genomes", fill = "Health Status") +
  theme_minimal() +  
  theme(
    legend.position = "bottom",
    legend.text = element_text(size=14),
    legend.title = element_blank(),
    axis.line = element_line(linewidth=0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size = 16, margin = margin(r = 10))
  )

resistance_plot <- ggplot(metadata, aes(x = factor(Resistance_score), fill = Health_Status)) +
  geom_histogram(stat="count", colour="grey", width=0.6, alpha=0.9) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +  
  labs(x = "Resistance score", y = "Number of genomes", fill = "Health Status") +
  theme_minimal() +  
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth=0.5),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size = 16, margin = margin(r = 10))
  )

# combine the two plots 
vir.amr.plot = ggarrange(virulence_plot, resistance_plot, ncol=2, common.legend = TRUE, legend = "bottom", labels=c("a", "b"), font.label = list(size=20))
ggsave(file="figures/vir-res_scores.pdf", height=5, width=10)

# wilcox tests
vir = wilcox.test(Virulence_score ~ Health_Status, data=metadata)
res = wilcox.test(Resistance_score ~ Health_Status, data=metadata)

vir.carriage = mean(metadata$Virulence_score[which(metadata$Health_Status == "Healthy")])
vir.disease = mean(metadata$Virulence_score[which(metadata$Health_Status == "Diseased")])

res.carriage = mean(metadata$Resistance_score[which(metadata$Health_Status == "Healthy")])
res.disease = mean(metadata$Resistance_score[which(metadata$Health_Status == "Diseased")])

# function to extract genes
extract_features = function(x, label) {
  subset.df = klebo[,which(grepl(x, colnames(klebo)))]
  subset.vec = as.vector(as.matrix(subset.df))
  subset.vec = subset.vec[which(subset.vec != "-")]
  subset.vec = unlist(strsplit(subset.vec, ";"))
  subset.fi = data.frame(table(subset.vec))
  subset.fi$Label = label
  return(subset.fi)
}
amr.acqu = extract_features("acquired", "AMR_genes")
amr.muta = extract_features("mutations", "AMR_mutations")

# plot frequency
amr.fi = rbind(amr.acqu, amr.muta)
amr.fi = amr.fi[order(amr.fi$Freq, decreasing=TRUE)[1:50],]
lolli = ggplot(amr.fi, aes(x=reorder(subset.vec, -Freq), y=Freq, fill=Label)) +
  geom_linerange(aes(ymin=0, ymax=Freq), colour="grey", linewidth=1) +
  geom_point(size=4, shape=21, colour="darkgrey") +
  theme_classic() +
  ylab("Number of genomes") +
  scale_fill_manual(values=c("darkolivegreen", "darkolivegreen2"), labels=c("Genes", "Mutations"), name="AMR determinant") +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) +
  theme(axis.title.y = element_text(size=16)) +
  theme(axis.text.y = element_text(size=14))

# combine final plot
ggarrange(vir.amr.plot, lolli, ncol=1, labels=c("", "c"), font.label = list(size=20))
ggsave("figures/vir-res_combined.pdf", height=9, width=12)

