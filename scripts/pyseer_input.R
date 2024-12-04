# load libraries
library(panstripe)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
pa = read.delim("panaroo/gene_presence_absence.Rtab")
rownames(pa) = pa$Gene
pa = pa[,-1]
metadata = read.delim("metadata/Metadata_14112024.tsv")
rownames(metadata) = metadata$Genome

# metadata filter
metadata = metadata[which(!is.na(metadata$Country) & !is.na(metadata$Health_Status)),]
metadata = metadata[which(metadata$Disease_Name %in% c("Infection", "Diarrhoea", "Healthy") & !is.na(metadata$Country)),]
pa = pa[,rownames(metadata)]
genes.keep = names(which(rowSums(pa) < ncol(pa)*0.9 & rowSums(pa) > ncol(pa)*0.01))
pa.fi = pa[genes.keep,]
pa.fi = cbind(rownames(pa.fi), pa.fi)
colnames(pa.fi)[1] = "Gene"

# save new table
write.table(pa.fi, file="gwas/test/pyseer-inf_input.tsv", quote=FALSE, row.names=FALSE, sep="\t")

