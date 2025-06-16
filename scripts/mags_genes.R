# load libraries
library(readxl)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
pa.raw = read.delim("panaroo/gene_presence_absence.Rtab", fileEncoding="UTF-8")
rownames(pa.raw) = pa.raw$Gene
pa.raw = pa.raw[,-1]
metadata = read.delim("metadata/Metadata_15022025.tsv")
rownames(metadata) = metadata$Genome

eggnog = read.delim("panaroo/eggnog_v2-1-3.tsv", check.names=FALSE)
deepvf = as.data.frame(read_xlsx("panaroo/mags-only_deepVF.xlsx"))
deeparg = read.delim("panaroo/deeparg.mapping.ARG", check.names = FALSE)
amrfinder = read.delim("panaroo/amrfinder_results.tsv")
arg.novel = deeparg$read_id[which(!deeparg$read_id %in% amrfinder$Contig.id)]

# prepare metadata files
metadata.gtype = metadata[,c("Genome", "Genome_Type")]
metadata.gtype$Genome_Type = ifelse(metadata.gtype$Genome_Type == "Isolate", 1, 0)

# COG dictionary
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

# find genes unique to MAGs
pa.raw = pa.raw[which(rowSums(pa.raw) > 0.01*ncol(pa.raw)),]
isolates = rownames(metadata.gtype)[which(metadata.gtype$Genome_Type == 1)]
mags = rownames(metadata.gtype)[which(metadata.gtype$Genome_Type == 0)]
isolat.genes = names(which(rowSums(pa.raw[,isolates]) > 0))
mags.genes = names(which(rowSums(pa.raw[,mags]) > 0))
mags.novel = mags.genes[which(!mags.genes %in% isolat.genes)]

# prepare final gene table
gene.df = data.frame(matrix(nrow=length(mags.novel),ncol=10))
colnames(gene.df) = c("Gene", "COG", "COG_category", "DeepARG", "DeepVF", "Healthy", "Diseased", "Unknown", "Total", "eggNOG_Annotation")
gene.df$Gene = mags.novel
gene.df$DeepARG = deeparg[match(gene.df$Gene, deeparg$read_id),"probability"]
gene.df$DeepVF = deepvf[match(gene.df$Gene, deepvf$Gene),"VF_score"]
gene.df$COG = eggnog[match(gene.df$Gene, eggnog[,1]), "COG_category"]
gene.df$COG = ifelse((gene.df$COG == "-" | is.na(gene.df$COG)), "S", gene.df$COG)
gene.df$COG_category = COG_descriptions[gene.df$COG]
gene.df$COG_category = ifelse(is.na(gene.df$COG_category), "Multiple functions", gene.df$COG_category)
gene.df$eggNOG_Annotation = eggnog[match(gene.df$Gene, eggnog[,1]), "Description"]
gene.df$eggNOG_Annotation = ifelse(gene.df$eggNOG_Annotation == "-", NA, gene.df$eggNOG_Annotation)

# count health/disease/all
for (gene in 1:nrow(gene.df)) {
  gene_name = gene.df$Gene[gene]
  select.genomes = colnames(pa.raw)[which(pa.raw[gene_name,] > 0)]
  gene.df$Total[gene] = length(select.genomes)
  gene.df$Healthy[gene] = length(which(metadata[select.genomes,"Health_Status"] == "Healthy"))
  gene.df$Diseased[gene] = length(which(metadata[select.genomes,"Health_Status"] == "Diseased"))
  gene.df$Unknown[gene] =  gene.df$Total[gene]-gene.df$Healthy[gene]- gene.df$Diseased[gene]
}

# plot
plot.df = gene.df[1:25,]
plot.df$COG_description = COG_descriptions[plot.df$COG]
plot.df$VF_class = ifelse(plot.df$DeepVF > 0.8, "Yes", "No")
bar.plot = ggplot(plot.df, aes(x=reorder(Gene, -Total), y=Total, fill=COG_description)) +
  geom_bar(stat="identity", width=0.6) +
  theme_classic() +
  ylab("Number of MAGs") +
  xlab("") +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, vjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.position = "bottom", legend.text = element_text(size=10), legend.title = element_text(size=14)) +
  scale_fill_manual(values=colorRampPalette(rev(brewer.pal(10, "Set3")))(length(unique(plot.df$COG_description))), name="COG category")

vf.plot = ggplot(plot.df, aes(x=reorder(Gene, -Total), y=DeepVF, fill=VF_class)) +
  geom_linerange(aes(ymin=0, ymax=DeepVF), colour="grey", linewidth=1) +
  geom_point(size=4, shape=21, colour="darkgrey") +
  geom_hline(yintercept=0.8, linetype = "dashed") +
  theme_classic() +
  xlab("") +
  ylab("DeepVF score") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.position = "top", legend.text = element_text(size=12), legend.title = element_text(size=14)) +
  scale_fill_manual(name="Putative virulence factor", values=c("steelblue", "tomato"))
ggarrange(vf.plot, bar.plot, nrow=2, heights=c(0.7,1), align="v",labels=c("a", "b"), font.label = list(size=18))
ggsave(file="figures/mag_genes.pdf", height=8, width=12)

# save genes
write.table(gene.df, file="tables/SupplementaryTable2.tsv", row.names=FALSE, quote=FALSE, sep="\t")