# load libraries
library(dendextend)
library(ape)
library(phytools)
library(ggtree)
library(phangorn)
library(TreeTools)
library(cluster)
library(fossil)
library(vegan)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/itol/")
panaroo = read.tree("species/panaroo.nwk")
gubbins = read.tree("species/snippy.tre")

# unique tip labels in each tree
unique_to_panaroo = setdiff(panaroo$tip.label, gubbins$tip.label)
unique_to_gubbins = setdiff(gubbins$tip.label, panaroo$tip.label)

cat("Number of tips in Panaroo tree:", length(panaroo$tip.label), "\n")
cat("Number of tips in Gubbins tree:", length(gubbins$tip.label), "\n")

# common tip labels in both trees
common_labels = intersect(panaroo$tip.label, gubbins$tip.label)
cat("Number of common tips:", length(common_labels), "\n")

# trees are ultrametric and midpoint-rooted
panaroo_ultra = chronos(panaroo)
panaroo_rooted = midpoint.root(panaroo_ultra)

gubbins_ultra = chronos(gubbins)
gubbins_rooted = midpoint.root(gubbins_ultra)

# convert to dendrograms
panaroo_dend = as.dendrogram(panaroo_rooted)
gubbins_dend = as.dendrogram(gubbins_rooted)

# reorder the dendrograms to have the same label order
labels_order = sort(labels(panaroo_dend))
dend1 = dendextend::rotate(panaroo_dend, labels_order)
dend2 = dendextend::rotate(gubbins_dend, labels_order)

# tanglegram
pdf(file="../figures/tanglegram_gene-vs-snp.pdf", height=5, width=10)
tanglegram(dend1, dend2, color_lines = "skyblue3", main_left = "Core tree (gene-based)", main_right = "Core tree (SNP-based)", match_order_by_labels = TRUE, 
           common_subtrees_color_lines = FALSE, highlight_branches_lwd = FALSE, highlight_distinct_edges = FALSE, 
           lab.cex = NULL, cex.axis = 1.5, margin_outer = 4, lwd=0.5, margin_inner = 0.5)
dev.off()

# get distances and perform mantel test
panaroo_cophenetic = cophenetic(panaroo)
gubbins_cophenetic = cophenetic(gubbins)
gubbins_cophenetic_ordered = gubbins_cophenetic[rownames(panaroo_cophenetic), colnames(panaroo_cophenetic)]
mantel_result = mantel(panaroo_cophenetic, gubbins_cophenetic_ordered, method="pearson", permutations=999)
