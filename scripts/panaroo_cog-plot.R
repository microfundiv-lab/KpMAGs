# load libraries
library(VennDiagram)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(viridis)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
core_genes = read.table("panaroo/core_genes_list.txt", header = FALSE, stringsAsFactors = FALSE)
eggnog_data = read.delim("panaroo/eggnog_v2-1-3.tsv", header = TRUE, sep = "\t")
colnames(eggnog_data)[1] = "query"

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

# split cogs
cog_data = eggnog_data[,c("query", "COG_category")]
cog_data$COG_category = gsub("-", "S", cog_data$COG_category)
expanded_cog <- cog_data %>%
  mutate(COG_category = strsplit(as.character(COG_category), "")) %>%
  unnest(COG_category) 
expanded_cog$COG_name = COG_descriptions[expanded_cog$COG_category]

# core genes data
pangenome = 20882
colnames(core_genes)[1] = "query"
core_total = nrow(core_genes)
acc_total = pangenome-core_total

# split data into core and accessory genomes
core_genome = expanded_cog[expanded_cog$query %in% core_genes$query, ]
accessory_genome = expanded_cog[!expanded_cog$query %in% core_genes$query, ]
core_gene_count = nrow(core_genome)
accessory_gene_count = nrow(accessory_genome)

# group core and accessory genomes by COG category
core_grouped = core_genome %>%
  group_by(COG_name) %>%
  summarize(Core_Count = n())

accessory_grouped = accessory_genome %>%
  group_by(COG_name) %>%
  summarize(Accessory_Count = n())

combined_grouped = full_join(core_grouped, accessory_grouped, by = "COG_name")

# replace NA values with 0 in Core_Count and Accessory_Count columns
combined_grouped <- combined_grouped %>%
  mutate(across(c(Core_Count, Accessory_Count), ~ replace_na(.x, 0)))

# normalize the gene counts by total gene counts for core and accessory genomes
combined_grouped <- combined_grouped %>%
  mutate(Core_Percentage = (Core_Count / core_gene_count) * 100,
         Accessory_Percentage = (Accessory_Count / accessory_gene_count) * 100)

# reshape the data for plotting
combined_grouped_long <- combined_grouped %>%
  pivot_longer(cols = c(Core_Percentage, Accessory_Percentage),
               names_to = "Genome_Type",
               values_to = "Gene_Percentage")


# calculate percentages per type
combined_grouped_long_bar = combined_grouped %>%
  pivot_longer(cols = c(Core_Count, Accessory_Count),
               names_to = "Genome_Type",
               values_to = "Gene_Count") %>%
  group_by(Genome_Type) %>%
  mutate(Gene_Percentage = (Gene_Count / sum(Gene_Count)) * 100)  # Convert to percentage

# find sig categories
counts.df = combined_grouped_long_bar
result.df = data.frame(variable=character(), core = numeric(), acc = numeric(), pvalue=numeric(), result=numeric())
var.list = unique(combined_grouped_long_bar$COG_name)
for (variable in var.list){
  counts.df$variable.group = ifelse(counts.df$COG_name == variable, variable, "all_other")
  cont.table = with(counts.df, xtabs(Gene_Count ~ Genome_Type + variable.group)) # create contingency table
  test.result = fisher.test(cont.table)
  core_prop = cont.table["Core_Count",variable]/(cont.table["Core_Count","all_other"]+cont.table["Core_Count",variable])*100
  acc_prop = cont.table["Accessory_Count",variable]/(cont.table["Accessory_Count","all_other"]+cont.table["Accessory_Count",variable])*100
  result = core_prop - acc_prop
  result.df = rbind(result.df, data.frame(variable = variable, core = core_prop, acc = acc_prop, pvalue = test.result$p.value, result = result))
  counts.df$variable.group = NULL # clear for next cycle
}
result.df$FDR = p.adjust(result.df$pvalue, method="fdr")
result.df = result.df[which(result.df$FDR < 0.05),]
result.df$Classification = ifelse(result.df$result > 0, "Core", "Accessory")

# plot sig categories
sig.cats = combined_grouped_long_bar[which(combined_grouped_long_bar$COG_name %in% result.df$variable),]
sig.cats = sig.cats[which(sig.cats$COG_name != "No functional prediction"),]

order.cats = result.df[order(result.df$result, decreasing=TRUE),"variable"]
order.cats = order.cats[which(order.cats != "No functional prediction")]

cog.plot = ggplot(sig.cats, aes(x = COG_name, y = Gene_Percentage, fill = Genome_Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "grey80", width = 0.7) + 
  coord_flip() + 
  labs(x = "COG category", y = "Proportion of genes (%)", fill = "Genome type") +
  scale_fill_manual(values = c("Core_Count" = "#377eb8", "Accessory_Count" = "#e41a1c"),
                    labels = c("Core_Count" = "Core genome", "Accessory_Count" = "Accessory genome")) +  
  scale_x_discrete(limits=order.cats) +
  theme_minimal() +  
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.title = element_blank(),  
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5),  
    panel.grid.minor = element_blank(),  
    panel.grid.major.y = element_blank(), 
    axis.line = element_line(linewidth=0.5),
    axis.ticks = element_blank(),  
    axis.title.x = element_text(size = 16, margin = margin(t = 15)),  
    axis.title.y = element_blank()   
  )

# calculate prop of unknown genes
core_match = length(unique(core_genome$query))
core_known = length(unique(core_genome[which(core_genome$COG_category != "S"),]$query))
acc_match = length(unique(accessory_genome$query))
acc_known = length(unique(accessory_genome[which(accessory_genome$COG_category != "S"),]$query))
core_nomatch = (core_total-core_match)/core_total*100
core_unknown = (core_total-core_known)/core_total*100
acc_nomatch = (acc_total-acc_match)/acc_total*100
acc_unknown = (acc_total-acc_known)/acc_total*100
