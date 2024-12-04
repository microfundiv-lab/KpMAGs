# load libraries
library(ggplot2)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
diseased.all = read.delim("gwas/test/pyseer-all_results-genes.tsv")
diseased.all$FDR = p.adjust(diseased.all$lrt.pvalue, method="fdr")
diseased.all = diseased.all[which(diseased.all$FDR < 0.05),]
diseased.inf = read.delim("gwas/test/pyseer-inf_results-genes.tsv")
diseased.inf$FDR = p.adjust(diseased.inf$lrt.pvalue, method="fdr")
diseased.inf = diseased.inf[which(diseased.inf$FDR < 0.05),]

# check overlap
diseased.all.pos = diseased.all$variant[which(diseased.all$beta > 0)]
diseased.inf.pos = diseased.inf$variant[which(diseased.inf$beta > 0)]
overlap.pos = intersect(diseased.all.pos, diseased.inf.pos)

diseased.all.neg = diseased.all$variant[which(diseased.all$beta < 0)]
diseased.inf.neg = diseased.inf$variant[which(diseased.inf$beta < 0)]
overlap.neg = intersect(diseased.all.neg, diseased.inf.neg)

# save output
write.table(overlap.pos, file="gwas/test/overlap-pos.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(overlap.neg, file="gwas/test/overlap-neg.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

