# load libraries
library(ape)
library(phytools)
library(picante)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/")
metadata = read.delim("metadata/Metadata_14112024.tsv")
rownames(metadata) = metadata$Genome
tree = read.tree("itol/species/panaroo.nwk")
tree = midpoint.root(tree)

#subsetting isolates & mags
metadata = metadata[which(metadata$Health_Status == "Healthy"),]
isolates = rownames(metadata)[which(metadata$Genome_Type == "Isolate")]
mags = rownames(metadata)[which(metadata$Genome_Type == "MAG")]

#subsetting the tree for isolates, MAGs, and combined
tree_isolates = drop.tip(tree, setdiff(tree$tip.label, isolates))
tree_mags = drop.tip(tree, setdiff(tree$tip.label, mags))
tree_combined = drop.tip(tree, setdiff(tree$tip.label, c(isolates, mags)))

#calculating Phylogenetic Diversity (PD) directly from the tree
pd_isolates_direct = sum(tree_isolates$edge.length)
pd_mags_direct = sum(tree_mags$edge.length)
pd_combined_direct = sum(tree_combined$edge.length)

#community matrix for isolates
comm_isolates = as.data.frame(matrix(1, nrow = 1, ncol = length(isolates)))
colnames(comm_isolates) = isolates
rownames(comm_isolates) = "Isolates"

#community matrix for MAGs
comm_mags = as.data.frame(matrix(1, nrow = 1, ncol = length(mags)))
colnames(comm_mags) = mags
rownames(comm_mags) = "MAGs"

#community matrix for the combined sample (Isolates + MAGs)
comm_combined = as.data.frame(matrix(1, nrow = 1, ncol = length(c(isolates, mags))))
colnames(comm_combined) = c(isolates, mags)
rownames(comm_combined) = "Combined"

faith_pd_isolates = pd(comm_isolates, tree_isolates)$PD
faith_pd_mags = pd(comm_mags, tree_mags)$PD
faith_pd_combined = pd(comm_combined, tree_combined)$PD

fc_pd = faith_pd_combined/faith_pd_isolates
fc_increas = (faith_pd_combined-faith_pd_isolates)/faith_pd_isolates*100
