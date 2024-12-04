# load libraries
library(ggplot2)
library(dplyr)
library(reshape2)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/")
pyseer_mags = read.delim("gwas/mags_vs_isolates//pyseer_results-genes.tsv", header = TRUE, sep = "\t")
pyseer_health = read.delim("gwas/diseased-inf/pyseer_results-genes.tsv", header = TRUE, sep = "\t")
eggnog = read.delim("panaroo/eggnog_v2-1-3.tsv", header = TRUE, sep = "\t")
colnames(eggnog)[1] = "variant"

# adjust pvalues
adjust_pvalues_and_merge = function(pyseer_df, eggnog_df) {
  pyseer_df$p.adjust = p.adjust(pyseer_df$lrt.pvalue, method = "fdr")
  merged_df = merge(pyseer_df, eggnog_df[, c("variant", "COG_category", "Description")], by = "variant", all.x = TRUE)
  return(merged_df)
}
pyseer_health = adjust_pvalues_and_merge(pyseer_health, eggnog)
pyseer_mags = adjust_pvalues_and_merge(pyseer_mags, eggnog)

# exclude mag genes
to_exclude = pyseer_mags[which(pyseer_mags$p.adjust < 0.05),"variant"]
pyseer_health = pyseer_health[which(!pyseer_health$variant %in% to_exclude),]

# sig cogs
process_data = function(data) {
  data$beta = as.numeric(as.character(data$beta))
  significant_results = data[data$p.adjust < 0.05, ]
  
  #split dataset based on the sign of beta
  positive_beta = significant_results[significant_results$beta > 0, ]
  negative_beta = significant_results[significant_results$beta < 0, ]
  
  #count different cog values in both subsets
  positive_counts = table(positive_beta$COG_category)
  negative_counts = table(negative_beta$COG_category)
  
  #tables to data frames
  positive_df = as.data.frame(positive_counts)
  negative_df = as.data.frame(negative_counts)
  
  names(positive_df) = c("COG", "Positive")
  names(negative_df) = c("COG", "Negative")
  
  #merge
  comparison_table = merge(positive_df, negative_df, by = "COG", all = TRUE)
  
  #replace NAs with zeros
  comparison_table$Positive[is.na(comparison_table$Positive)] = 0
  comparison_table$Negative[is.na(comparison_table$Negative)] = 0
  
  return(comparison_table)
}
COG_pyseer_health = process_data(pyseer_health)

# remove '-' and S
COG_pyseer_health$COG = gsub("-", "S", COG_pyseer_health$COG)

split_COGs_and_add_counts = function(df) {

  df$COG = as.character(df$COG)
  unique_cogs = unique(unlist(strsplit(unique(df$COG), "")))
  
  #initialize a df to hold all unique COGs with zero counts
  all_cogs_df = data.frame(COG = unique_cogs, Positive = numeric(length(unique_cogs)), Negative = numeric(length(unique_cogs)), stringsAsFactors = FALSE)
  all_cogs_df$Positive = 0
  all_cogs_df$Negative = 0
  
  #update counts for single-letter COGs 
  single_letter_rows = nchar(df$COG) == 1
  for (cog in df$COG[single_letter_rows]) {
    all_cogs_df[all_cogs_df$COG == cog, "Positive"] = sum(df$Positive[df$COG == cog])
    all_cogs_df[all_cogs_df$COG == cog, "Negative"] = sum(df$Negative[df$COG == cog])
  }
  
  #process combined COG entries
  combined_letter_rows = nchar(df$COG) > 1
  for (i in which(combined_letter_rows)) {
    row = df[i, ]
    letters = strsplit(row$COG, "")[[1]]
    
    #add counts to each constituent letter
    for (letter in letters) {
      all_cogs_df[all_cogs_df$COG == letter, "Positive"] = all_cogs_df[all_cogs_df$COG == letter, "Positive"] + row$Positive
      all_cogs_df[all_cogs_df$COG == letter, "Negative"] = all_cogs_df[all_cogs_df$COG == letter, "Negative"] + row$Negative
    }
  }
  return(all_cogs_df)
}
COG_pyseer_health = split_COGs_and_add_counts(COG_pyseer_health)

# fisher test
COG_pyseer_long = reshape2::melt(COG_pyseer_health)
result.df = data.frame(variable=character(), core = numeric(), acc = numeric(), pvalue=numeric(), result=numeric())
var.list = unique(COG_pyseer_long$COG)
for (cog.code in var.list){
  COG_pyseer_long$group = ifelse(COG_pyseer_long$COG == cog.code, cog.code, "all_other")
  cont.table = with(COG_pyseer_long, xtabs(value ~ variable + group)) # create contingency table
  test.result = fisher.test(cont.table)
  pos_prop = cont.table["Positive",cog.code]/(cont.table["Positive","all_other"]+cont.table["Positive",cog.code])*100
  neg_prop = cont.table["Negative",cog.code]/(cont.table["Negative","all_other"]+cont.table["Negative",cog.code])*100
  result = pos_prop - neg_prop
  result.df = rbind(result.df, data.frame(variable = cog.code, pos = pos_prop, neg = neg_prop, pvalue = test.result$p.value, result = result))
  COG_pyseer_long$group = NULL # clear for next cycle
}
result.df$FDR = p.adjust(result.df$pvalue, method="fdr")
result.df = result.df[which(result.df$FDR < 0.2),]
result.df$Classification = ifelse(result.df$result > 0, "Diseased", "Healthy")

# save file
write.csv(COG_pyseer_health, "gwas/diseased-inf/pyseer_results-cogs_nomags.csv", row.names = FALSE)
