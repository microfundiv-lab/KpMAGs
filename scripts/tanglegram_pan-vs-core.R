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
library(panstripe)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
panaroo <- read.tree("itol/species/panaroo.nwk")
panaroo.ultra <- chronos(panaroo)
panaroo.rooted <- midpoint.root(panaroo.ultra)
pan_dist = read_rtab("panaroo/gene_presence_absence.Rtab")
select.genes = names(which(colSums(pan_dist) > nrow(pan_dist)*0.01 & colSums(pan_dist) < nrow(pan_dist)*0.9))
pangenome = vegdist(pan_dist[,select.genes], method="jaccard")
pangenome.hc = hclust(pangenome)

# convert to dendrograms
panaroo_dend = as.dendrogram(panaroo.rooted)
pangenome_dend = as.dendrogram(pangenome.hc)

# reorder the dendrograms to have the same label order
labels_order <- sort(labels(panaroo_dend))
dend1 <- dendextend::rotate(panaroo_dend, labels_order)
dend2 <- dendextend::rotate(pangenome_dend, labels_order)

# tanglegram
pdf(file="figures/tanglegram_pan-vs-core.pdf", height=5, width=10)
tanglegram(dend1, dend2, color_lines = "skyblue3", main_left = "Core tree (gene-based)", main_right = "Pan-genome distances", match_order_by_labels = TRUE, 
           common_subtrees_color_lines = FALSE, highlight_branches_lwd = FALSE, highlight_distinct_edges = FALSE, 
           lab.cex = NULL, cex.axis = 1.5, margin_outer = 4, lwd=0.5, margin_inner = 0.5)
dev.off()

# get distances and perform mantel test
panaroo_cophenetic <- cophenetic(panaroo)
pangenome_dist <- as.matrix(pangenome)
pangenome_ordered <- pangenome_dist[rownames(panaroo_cophenetic), colnames(panaroo_cophenetic)]
mantel_result <- mantel(panaroo_cophenetic, pangenome_ordered, method="pearson", permutations=999)
