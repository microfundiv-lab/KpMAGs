# load libraries
library(ape)
library(phytools)
library(picante)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/")
metadata = read.delim("tables/SupplementaryTable1.tsv")
rownames(metadata) = metadata$Genome
tree = read.tree("itol/species/panaroo.nwk")
tree = midpoint.root(tree)

# subsetting isolates & mags
isolate.countries = unique(metadata[which(metadata$Genome.type == "Isolate"),"Country"])
mags.countries = unique(metadata[which(metadata$Genome.type == "MAG" & metadata$Mash > 0.005),'Country'])
common.countries = intersect(isolate.countries, mags.countries)

# subset files by country
get_pd_country = function(country) {
  metadata.cont = metadata[which(metadata$Country == country),]
  tree.cont = drop.tip(tree, setdiff(tree$tip.label, rownames(metadata.cont)))
  isolates = rownames(metadata.cont)[which(metadata.cont$Genome.type == "Isolate")]
  mags = rownames(metadata.cont)[which(metadata.cont$Genome.type == "MAG" & metadata.cont$Mash > 0.005)]
  
  # subsetting the tree for isolates, MAGs, and novel
  tree_isolates = drop.tip(tree.cont, setdiff(tree$tip.label, isolates))
  tree_mags = drop.tip(tree.cont, setdiff(tree$tip.label, mags))
  
  # calculating Phylogenetic Diversity (PD)
  pd_isolates = sum(tree_isolates$edge.length)
  pd_mags = sum(tree_mags$edge.length)
  pd_all = sum(tree.cont$edge.length)
  
  # calculate improvement
  mags.increase = (pd_all-pd_isolates)/pd_isolates*100
  mags.fc = pd_all/pd_isolates
  return(c(mags.increase, mags.fc))
}

# group data
count.df = data.frame(matrix(nrow=length(common.countries), ncol=2, NA))
rownames(count.df) = common.countries
colnames(count.df) = c("Increase%", "FC")
for (c in common.countries){
  count.df[c,] = get_pd_country(c)
}
count.df = count.df[order(count.df$FC, decreasing=TRUE),]
count.df$Country = rownames(count.df)

# plot
phylo = ggplot(count.df, aes(x=reorder(Country, FC), y=FC, fill=Country)) +
  geom_linerange(aes(ymin=0, ymax=FC), colour="grey", linewidth=1) +
  geom_point(size=4, shape=21, colour="darkgrey") +
  coord_flip() +
  theme_classic() +
  xlab("") +
  ylab("PD increase\n(Fold change)") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.position = "none", legend.text = element_text(size=12), legend.title = element_text(size=14)) +
  scale_fill_manual(name="Country", values=colors.dict)
