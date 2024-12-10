# load libraries
library(panstripe)
library(ape)
library(patchwork)
library(phytools)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(gplots)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
pa = read_rtab("panaroo/gene_presence_absence.Rtab")
tree = read.tree("itol/species/panaroo.nwk")
tree = midpoint.root(tree)
metadata = read.delim("metadata/Metadata_09122024.tsv")
rownames(metadata) = metadata$Genome
gwas.mags = read.delim("gwas/gtype_results-genes.tsv")
gwas.mags$FDR = p.adjust(gwas.mags$lrt.pvalue, method="fdr")
gwas.mags = gwas.mags[which(gwas.mags$FDR < 0.05),]
gwas = read.delim("gwas/diseased-inf_results-genes.tsv")
gwas$FDR = p.adjust(gwas$lrt.pvalue, method="fdr")
gwas = gwas[which(gwas$FDR < 0.05 & !gwas$variant %in% gwas.mags$variant),]
gwas.pos = gwas[which(gwas$beta > 0),]
gwas.neg = gwas[which(gwas$beta < 0),]

# metadata filter
metadata = metadata[which(metadata$Disease_Name %in% c("Infection", "Diarrhoea", "Healthy") & !is.na(metadata$Country)),]
toremove = rownames(pa)[which(!rownames(pa) %in% rownames(metadata))]
pa = pa[rownames(metadata),]
tree = drop.tip(tree, toremove)

# create itol dataset with country
itol.df = metadata
itol.country = itol.df[,c("Genome", "Country")]
itol.country = itol.country[!is.na(itol.country$Country),]
top10 = names(sort(table(itol.country$Country), decreasing=TRUE)[1:9])
itol.country$Country = ifelse(itol.country$Country %in% top10, itol.country$Country, "Other")
colors = c(brewer.pal(9,"Set1"), col2hex("grey50"))
names(colors) = c(top10, "Other")
itol.country$Color = as.vector(colors[itol.country$Country])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,75\nSTRIP_WIDTH,150\nDATASET_LABEL,Country\nSIZE_FACTOR,7\nLABEL_SHIFT,20\nSHOW_LABELS,1\nDATA\n", file="itol/gwas-inf/itol_country.txt")
write.table(itol.country[,c("Genome", "Color", "Country")], file="itol/gwas-inf/itol_country.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with ST
itol.st = itol.df[,c("Genome", "ST")]
itol.st = itol.st[!is.na(itol.st$ST),]
top = names(sort(table(itol.st$ST), decreasing=TRUE)[1:12])
itol.st$ST = ifelse(itol.st$ST %in% top, itol.st$ST, "Other")
colors = c(brewer.pal(12,"Set3"), col2hex("grey50"))
names(colors) = c(top, "Other")
itol.st$Color = as.vector(colors[itol.st$ST])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,75\nSTRIP_WIDTH,150\nDATASET_LABEL,ST\nSIZE_FACTOR,7\nLABEL_SHIFT,20\nSHOW_LABELS,1\nDATA\n", file="itol/gwas-inf/itol_st.txt")
write.table(itol.st[,c("Genome", "Color", "ST")], file="itol/gwas-inf/itol_st.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with health state
itol.state = itol.df[,c("Genome", "Health_Status")]
itol.state = itol.state[!is.na(itol.state$Health_Status),]
colors = c("#0072B2", "#D55E00")
names(colors) = c("Healthy", "Diseased")
itol.state$Color = as.vector(colors[itol.state$Health_Status])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,75\nSTRIP_WIDTH,150\nDATASET_LABEL,Health status\nSIZE_FACTOR,7\nLABEL_SHIFT,20\nSHOW_LABELS,1\nDATA\n", file="itol/gwas-inf/itol_status.txt")
write.table(itol.state[,c("Genome", "Color", "Health_Status")], file="itol/gwas-inf/itol_status.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with source
itol.source = itol.df[,c("Genome", "Source")]
itol.source = itol.source[!is.na(itol.source$Source),]
colors = c("#ffd966", "#cfe2f3")
names(colors) = c("Community", "Hospital")
itol.source$Color = as.vector(colors[itol.source$Source])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,75\nSTRIP_WIDTH,150\nDATASET_LABEL,Source\nSIZE_FACTOR,7\nLABEL_SHIFT,20\nSHOW_LABELS,1\nDATA\n", file="itol/gwas-inf/itol_source.txt")
write.table(itol.source[,c("Genome", "Color", "Source")], file="itol/gwas-inf/itol_source.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with genome type
itol.gtype = itol.df[,c("Genome", "Genome_Type")]
itol.gtype$Color = ifelse(itol.gtype$Genome_Type == "Isolate", "range,#b4c8ff,Isolate", "range,#d9ffc1,MAG")
cat("TREE_COLORS\nSEPARATOR COMMA\nLEGEND_TITLE,Genome type\nDATA\n", file="itol/gwas-inf/itol_gtype.txt")
write.table(itol.gtype[,-2], file="itol/gwas-inf/itol_gtype.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with gwas results
itol.pos = pa[,gwas.pos$variant]
itol.neg = pa[,gwas.neg$variant]
itol.neg[itol.neg == 1] = -1
itol.gwas = cbind(itol.pos, itol.neg)
itol.gwas = itol.gwas[,names(sort(colSums(abs(itol.gwas)), decreasing=TRUE))]
itol.gwas = cbind(rownames(itol.gwas), itol.gwas)
colnames(itol.gwas)[1] = "Genome"
cat(paste0("DATASET_HEATMAP\nSEPARATOR COMMA\nDATASET_LABEL,Effect size\nMARGIN,200\nSTRIP_WIDTH,20\nSIZE_FACTOR,0\nCOLOR_MIN,#0072B2\nCOLOR_MAX,#D55E00\nUSE_MID_COLOR,1\nCOLOR_MID,#ffffff\nUSER_MID_VALUE,0\nSHOW_LABELS,0\nFIELD_LABELS,", paste(colnames(itol.gwas[,-1]), collapse=","), "\nDATA\n"), file="itol/gwas-inf/itol_gwas.txt")
write.table(itol.gwas, file="itol/gwas-inf/itol_gwas.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# save tree
write.tree(tree, file="itol/gwas-inf/gwas-inf.nwk")