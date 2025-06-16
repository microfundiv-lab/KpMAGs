# load libraries
library(ggplot2)
library(phytools)
library(tidyr)
library(RColorBrewer)
library(gplots)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")

# load da data
itol.df = read.delim("metadata/Metadata_15022025.tsv", check.names=FALSE)
rownames(itol.df) = itol.df$Genome
novel.lineages = scan("mash/novel_lineages_mag.txt", what="")
itol.df = itol.df[novel.lineages,]

# create itol dataset with country
itol.country = itol.df[,c("Genome", "Country")]
itol.country = itol.country[!is.na(itol.country$Country),]
top10 = names(sort(table(itol.country$Country), decreasing=TRUE)[1:9])
itol.country$Country = ifelse(itol.country$Country %in% top10, itol.country$Country, "Other")
colors = c(brewer.pal(9,"Set1"), col2hex("grey50"))
names(colors) = c(top10, "Other")
itol.country$Color = as.vector(colors[itol.country$Country])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,40\nSTRIP_WIDTH,50\nDATASET_LABEL,Country\nSIZE_FACTOR,2\nLABEL_SHIFT,20\nLABEL_ALIGN_TO_TREE,0\nSHOW_LABELS,1\nDATA\n", file="itol/novel/itol_country.txt")
write.table(itol.country[,c("Genome", "Color", "Country")], file="itol/novel/itol_country.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with ST
itol.st = itol.df[,c("Genome", "ST")]
itol.st = itol.st[!is.na(itol.st$ST),]
top = names(sort(table(itol.st$ST), decreasing=TRUE)[1:12])
itol.st$ST = ifelse(itol.st$ST %in% top, itol.st$ST, "Other")
colors = c(brewer.pal(12,"Set3"), col2hex("grey50"))
names(colors) = c(top, "Other")
itol.st$Color = as.vector(colors[itol.st$ST])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,40\nSTRIP_WIDTH,50\nDATASET_LABEL,ST\nSIZE_FACTOR,2\nLABEL_SHIFT,20\nLABEL_ALIGN_TO_TREE,0\nSHOW_LABELS,1\nDATA\n", file="itol/novel/itol_st.txt")
write.table(itol.st[,c("Genome", "Color", "ST")], file="itol/novel/itol_st.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with health state
itol.state = itol.df[,c("Genome", "Health_Status", "Disease_Name")]
itol.state$Health_Status = ifelse(itol.state$Disease_Name %in% c("Infection", "Diarrhoea"), "Diseased (infection)", ifelse(!is.na(itol.state$Disease_Name), "Diseased (other)", itol.state$Disease_Name))
itol.state$Health_Status = ifelse(itol.state$Disease_Name == "Healthy", "Healthy", itol.state$Health_Status)
itol.state = itol.state[!is.na(itol.state$Health_Status),]
colors = c("#0072B2", "#cd5f08", "#ffa863")
names(colors) = c("Healthy", "Diseased (infection)", "Diseased (other)")
itol.state$Color = as.vector(colors[itol.state$Health_Status])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,40\nSTRIP_WIDTH,50\nDATASET_LABEL,Health status\nSIZE_FACTOR,2\nLABEL_SHIFT,20\nSHOW_LABELS,1\nLABEL_ALIGN_TO_TREE,0\nDATA\n", file="itol/novel/itol_status.txt")
write.table(itol.state[,c("Genome", "Color", "Health_Status")], file="itol/novel/itol_status.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)

# create itol dataset with source
itol.source = itol.df[,c("Genome", "Source")]
itol.source = itol.source[!is.na(itol.source$Source),]
colors = c("#ffd966", "#cfe2f3")
names(colors) = c("Community", "Hospital")
itol.source$Color = as.vector(colors[itol.source$Source])
cat("DATASET_COLORSTRIP\nSEPARATOR COMMA\nMARGIN,40\nSTRIP_WIDTH,50\nDATASET_LABEL,Source\nSIZE_FACTOR,2\nLABEL_SHIFT,20\nLABEL_ALIGN_TO_TREE,0\nSHOW_LABELS,1\nDATA\n", file="itol/novel/itol_source.txt")
write.table(itol.source[,c("Genome", "Color", "Source")], file="itol/novel/itol_source.txt", sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)