# load libraries
library(matrixStats)
library(panstripe)
library(ggplot2)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/")
pa = read_rtab("panaroo/gene_presence_absence.Rtab") 

# define function
panplots = function(data, curve = "pan", iterations = 5, thresh) {
 
  nr_rows = nrow(data);
  nr_iterations = iterations;
  
  #create empty matrix to store temp results
  temp = matrix(data=NA,nrow=nr_rows,ncol=nr_iterations)
  
  if(curve == "core") {
    
    ## compute core_genome_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      
      # random sampling of genomes
      for (i in 2: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) >= thresh*i))
      }
    }
  } 
  
  if(curve == "pan") {
    ## compute gene_cluster_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 1: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) > thresh*i))
      }
    }
  } 
  
  # summarize permutation results using "matrixStats" library
  summary = data.frame(genomes=c(1:nr_rows)) 
  summary$mean=rowMeans2(temp[,c(-1)])
  summary$sd=rowSds(temp[,c(-1)])
  summary$group=deparse(substitute(data))
  return(summary)
}

pan.core90 = panplots(pa, curve="core", iterations=5, thresh=0.90)
pan.core90$genes = "Core (90%)"

pan.core100 = panplots(pa, curve="core", iterations=5, thresh=1.00)
pan.core100$genes = "Core (100%)"

pan1 = panplots(pa, curve="pan", iterations=5, thresh=0.01)
pan1$genes = "Pan-genome (1%)"

pan.all = panplots(pa, curve="pan", iterations=5, thresh=0)
pan.all$genes = "Pan-genome (all)"

pan.fi = rbind(pan.core90, pan.core100, pan.all, pan1)
select.breaks = seq(1,nrow(pa), 50)
pan.fi = pan.fi[which(pan.fi$genomes %in% select.breaks),]
pan.fi$genes = factor(pan.fi$genes, levels=unique(pan.fi$genes))

# plot pangenome curves
pan.plot = ggplot(pan.fi, aes(x=genomes, y=mean, colour=genes, fill=genes)) +
  geom_point(alpha=0.8, size=1) +
  geom_line(linetype = "solid") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), linewidth=0.5, width=10) +
  ylab("Number of genes") +
  xlab("Number of genomes") +
  scale_colour_manual(values = c("steelblue", "steelblue1", "darkolivegreen", "darkolivegreen3"), name="Gene type") +
  scale_fill_manual(values = c("steelblue", "steelblue1", "darkolivegreen", "darkolivegreen3"), name="Gene type") +
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size=14)) +
  theme(axis.title = element_text(size=16)) +
  theme(axis.text = element_text(size=14))

# plot gene frequency
gene.freqs = data.frame(colSums(pa))
colnames(gene.freqs) = "Frequency"
freq.plot = ggplot(gene.freqs, aes(x=Frequency)) + 
  geom_histogram(fill="grey80", colour="grey40", linewidth=0.3) + 
  geom_vline(xintercept=nrow(pa)*0.9, linetype="dashed") +
  theme_classic() +
  ylab("Number of genes") +
  xlab("Number of genomes") +
  theme(axis.text = element_text(size=14)) +
  theme(axis.title = element_text(size=16))

# combine plots
source("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/scripts/alex/panaroo_cog-plot.R")
ggarrange(ggarrange(freq.plot, pan.plot, ncol=1, heights=c(1,1.5), labels=c("a", "b"), font.label = list(size=20)), cog.plot, ncol=2, labels=c("","c"), widths=c(1,1.5), font.label = list(size=20))
ggsave(file="figures/panaroo_combined.pdf", height=8, width=14)
