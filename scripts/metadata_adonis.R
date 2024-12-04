# load libraries
library(tidyverse)
library(ggplot2)
library(ggtree)
library(ape)
library(vegan)
library(data.table)
library(panstripe)

setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
metadata = read.delim("metadata/Metadata_14112024.tsv")
rownames(metadata) = metadata$Genome

# reformat health status
metadata$Health_Status = ifelse(metadata$Disease_Name %in% c("Infection", "Diarrhoea"), "Diseased (infection)", ifelse(!is.na(metadata$Disease_Name), "Diseased (other)", metadata$Disease_Name))
metadata$Health_Status = ifelse(metadata$Disease_Name == "Healthy", "Healthy", metadata$Health_Status)

gen_adonis = function(meta, data, type){  
  
  metadata.subs = metadata[!is.na(metadata[,meta]),] 
  if (type %in% c("Core tree (gene-based)", "Core tree (SNP-based)")) {
    phylo = read.tree(data)
    coph_dist = cophenetic(phylo)
    common_samples = intersect(rownames(coph_dist), metadata.subs$Genome)
    coph_dist = coph_dist[common_samples, common_samples]
    dist.fi= as.dist(coph_dist)
  } else if (type == "Pan-genome") {
    pan_dist = read_rtab(data)
    common_samples = intersect(rownames(pan_dist), metadata.subs$Genome)
    pan_dist = pan_dist[common_samples,]
    select.genes = names(which(colSums(pan_dist) > nrow(pan_dist)*0.01 & colSums(pan_dist) < nrow(pan_dist)*0.9))
    dist.fi = vegdist(pan_dist[,select.genes], method="jaccard")
  }
  metadata.subs = metadata.subs[common_samples,]

  # adonis
  adonis.test = adonis2(as.formula(paste0("dist.fi ~ ", meta)), data = metadata.subs)
  
  # save output
  df = data.frame(Data=type, Variable=meta, R2=adonis.test$R2[1], Pvalue=adonis.test$`Pr(>F)`[1], Df=adonis.test$Df[1], N=nrow(metadata.subs))
  return(df)
}

# generate effect per variable
variables = c("Continent", "Country", "Genome_Type", "Source", "Health_Status", "Disease_Name")
gen_adonis_multiple = function(data, type) {
  adonis.list = lapply(variables, function(x) {
    cat("Performing adonis for", x, "in", type, "...\n")
    adonis.data = gen_adonis(x, data, type)
    return(adonis.data)
  })
  adonis.combined = as.data.frame(rbindlist(adonis.list))
  adonis.combined$Type = type
  return(adonis.combined)
}

# run function for core and pan
panaroo_adonis = gen_adonis_multiple("itol/species/panaroo.nwk", "Core tree (gene-based)")
snippy_adonis = gen_adonis_multiple("itol/species/snippy.tre", "Core tree (SNP-based)")
pan_adonis = gen_adonis_multiple("panaroo/gene_presence_absence.Rtab", "Pan-genome")

# combine data
adonis_all = rbind(panaroo_adonis, snippy_adonis, pan_adonis)

# adjust R2
adonis_all$R2_adj = 1 - ((1 - adonis_all$R2) * (adonis_all$N - 1) / (adonis_all$N - adonis_all$Df - 1))
adonis_all$R2_adj = adonis_all$R2_adj*100
adonis_all$R2 = adonis_all$R2*100
adonis_all$Variable = gsub("_", " ", adonis_all$Variable)
adonis_all$Variable = gsub("Disease Name", "Health status\n(specific)", adonis_all$Variable)
adonis_all$Variable = gsub("Genome Type", "Genome type", adonis_all$Variable)
adonis_all$Variable = gsub("Health Status", "Health status\n(generic)", adonis_all$Variable)

# plot
adonis.plot = ggplot(adonis_all, aes(x = reorder(Variable,R2_adj), y = R2_adj, fill = Type)) +
  geom_bar(stat = "identity", position=position_dodge(width = 0.8), width = 0.8, color = "grey80", alpha=0.8) +  
  scale_fill_manual(values = c("darkolivegreen3", "darkolivegreen", "#377eb8")) +  
  ylab(expression("Effect size (R"^2~", %)"))+ 
  geom_text(aes(label = sprintf("%.0f%%", R2_adj)), position = position_dodge(width = 0.8), hjust = -0.1, size = 4.25, color = "black") +  
  coord_flip() +
  theme_classic() +  
  theme(
    axis.title.x = element_text(size = 14, margin = margin(r = 10)),  
    axis.title.y = element_blank(),  
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12), 
    plot.title = element_blank(),  
    axis.line = element_line(color = "black", linewidth=0.5),  
    axis.ticks = element_blank(),
    legend.box.margin = margin(0,0,0,-100),
    legend.position = "inside",
    legend.position.inside=c(0.9,0.3), 
    legend.box = "horizontal", 
    legend.text = element_text(size=12),
    legend.title = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
ggsave(file = "figures/metadata_adonis.pdf", height=8, width=4)
