# load libraries
library(panstripe)
library(ape)
library(patchwork)
library(phytools)
library(vegan)
library(ggplot2)

set.seed(1234)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
pa = read_rtab("panaroo/gene_presence_absence.Rtab")
tree = read.tree("itol/panaroo.nwk")
tree = midpoint.root(tree)
metadata = read.delim("metadata/Metadata_14112024.tsv")
rownames(metadata) = metadata$Genome

# metadata filter
#metadata = metadata[which(metadata$Disease_Name %in% c("Infection", "Diarrhoea", "Healthy") & !is.na(metadata$Country)),]
#toremove = rownames(pa)[which(!rownames(pa) %in% rownames(metadata))]
#pa = pa[rownames(metadata),]
#tree = drop.tip(tree, toremove)

# run panstripe
fit = panstripe(pa, tree, fit_method = "glm", family = "poisson")
fit$summary

# plot summaries
plot_pangenome_params(fit)
plot_pangenome_cumulative(fit)
plot_tsne(pa, category=metadata[rownames(pa),"Health_Status"])
plot_acc(pa)

# plot tree with accessory genes
plot_tree_pa(tree = tree, pa = pa, genes = sig.genes, label_genes = FALSE, cols = "black")

# plot gene frequency
gene.freqs = data.frame(colSums(pa))
colnames(gene.freqs) = "Frequency"
freq.plot = ggplot(gene.freqs, aes(x=Frequency)) + 
  geom_histogram(fill="grey80", colour="grey40", linewidth=0.3) + 
  geom_vline(xintercept=nrow(pa)*0.9, linetype="dashed") +
  theme_classic() +
  ylab("Number of genes") +
  xlab("Number of genomes") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

# pca
select.genes = names(which(colSums(pa) > nrow(pa)*0.05))
gene.dist = vegdist(pa[,select.genes], method="jaccard")
pcoa.data = pcoa(gene.dist)
pcoa.axes = data.frame(pcoa.data$vectors)
pc1_var = round(pcoa.data$values[1,2]*100,1)
pc2_var = round(pcoa.data$values[2,2]*100,1)
pca.df = data.frame(row.names=rownames(pcoa.axes), PC1=pcoa.axes[,1], PC2=pcoa.axes[,2],
                    Health_Status=metadata[rownames(pcoa.axes),"Health_Status"],
                    Genome_Type=metadata[rownames(pcoa.axes),"Genome_Type"],
                    Country=metadata[rownames(pcoa.axes),"Country"])

pca.df = pca.df[which(pca.df$Health_Status != "NA"),]

pca.plot = ggplot(pca.df, aes(x=PC1, y=PC2, colour=Health_Status)) + 
  geom_point(size=0.5, alpha=1) +
  theme_classic() +
  #stat_ellipse() +
  theme(panel.grid.minor = element_blank()) + 
  theme(legend.pos = "top") +
  ylab(paste("PC2"," (",pc2_var,"%)",sep="")) + 
  xlab(paste("PC1"," (",pc1_var,"%)",sep="")) +
  scale_colour_manual(values=c("tomato", "steelblue"), name="Health status") +
  theme(legend.position="top", legend.text=element_text(size=12),
        legend.title = element_text(size=12, face="bold")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.title.y = element_text(size=14)) + 
  theme(axis.text.y = element_text(size=12)) + 
  theme(axis.title.x = element_text(size=14)) + 
  theme(axis.text.x = element_text(size=12))

# adonis2
gene.df = as.data.frame(as.matrix(gene.dist))
gene.df = gene.df[rownames(metadata), rownames(metadata)]
perm.test = adonis2(gene.df ~ Country + Genome_Type + Health_Status, data=metadata)
