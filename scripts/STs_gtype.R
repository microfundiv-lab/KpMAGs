# load libraries
library(readr)
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster)
library(ape)
library(tidyverse)
library(glmmTMB)
library(patchwork)
library(ggvenn)
library(ggpubr)
library(grid)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/metadata")
cmseq = read.delim("CMseq.tsv")
combined_data = read.delim("Metadata_15022025.tsv")
combined_data = combined_data[which(combined_data$ST != "0"),]

# add LV column
split_matrix = do.call(rbind, strsplit(as.character(combined_data$ST), "-"))
combined_data$ST_1 = split_matrix[,1]
combined_data$ST_2 = split_matrix[,2]
combined_data$ST_2 = ifelse(grepl("LV", combined_data$ST_2), combined_data$ST_2, NA)

# subset countries
country_gtype_table = table(combined_data$Country, combined_data$Genome_Type)
country_gtype_matrix = as.data.frame.matrix(country_gtype_table)
select.countries = rownames(country_gtype_matrix)[which(country_gtype_matrix$MAG > 1 & country_gtype_matrix$Isolate > 1 )]
combined_data_filt = combined_data[which(combined_data$Country %in% select.countries),]
country_gtype_table = table(combined_data_filt$ST, combined_data_filt$Genome_Type)
country_gtype_matrix = as.data.frame.matrix(country_gtype_table)

# ST analysis
ST_gtype_table = table(combined_data$ST, combined_data$Genome_Type)
ST_gtype_matrix = as.data.frame.matrix(ST_gtype_table)
ST_gtype_df = as.data.frame(ST_gtype_table)

# check common STs
common.sts = names(which(rowSums(ST_gtype_matrix > 0) == 2))
common.df = combined_data[which(combined_data$ST %in% common.sts),]

result <- common.df %>%
  group_by(ST, Country) %>%
  filter(n_distinct(Genome_Type) > 1) %>%
  ungroup()

# extract top
top_isolate_STs = combined_data %>%
  filter(Genome_Type == "Isolate") %>%
  group_by(ST) %>%
  summarize(Count = n()) %>%
  arrange(desc(Count)) %>%
  top_n(8, Count)

top_MAG_STs = combined_data %>%
  filter(Genome_Type == "MAG") %>%
  group_by(ST) %>%
  summarize(Count = n()) %>%
  arrange(desc(Count)) %>%
  top_n(5, Count)

# plot for isolates
plot_isolate = ggplot(top_isolate_STs, aes(x = reorder(ST, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = '#6a92ff', width = 0.4, alpha=0.8, colour="darkgrey") +
  labs(x = "Sequence Type (ST)", y = "") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),  
    axis.title.x = element_text(size = 18), 
    axis.text.x = element_text(size = 14, angle=45, hjust=1, vjust=1),
    axis.text.y = element_text(size = 14),
    plot.margin=unit(c(1.5,1,1,1), 'cm'),
    plot.title = element_blank(),  
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dashed'),
    panel.grid.minor = element_blank() 
  ) 

# plot for MAGs
plot_MAG = ggplot(top_MAG_STs, aes(x = reorder(ST, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = '#90e25c', width = 0.4, alpha=0.6, colour="darkgrey") +
  labs(x = "", y = "") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),  
    axis.title.x = element_text(size = 18),  
    axis.text.x = element_text(size = 14, angle=45, hjust=1, vjust=1),  
    axis.text.y = element_text(size = 14),  
    plot.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dashed'),
    panel.grid.minor = element_blank() 
  )

combined_plot = ggarrange(plot_MAG, plot_isolate, align="v", ncol=1, heights=c(1,1.2))
ggsave(file="../figures/STs_gtype_bar.pdf", width=4, height=6)

# get values for venn diagram
gen_venn_list = function(x) {
  isolate.list = rownames(x)[which(x$Isolate > 0)]
  mag.list = rownames(x)[which(x$MAG > 0)]
  venn.list = list(Isolate = isolate.list, MAG = mag.list)
  return(venn.list)
}

venn.all = ggvenn(gen_venn_list(ST_gtype_matrix), fill_color = c("#6a92ff", "#90e25c"), stroke_size = 0, set_name_size = 0, auto_scale = TRUE, text_size=6, show_outside = "none", padding=0.15) + coord_flip()
venn.filt = ggvenn(gen_venn_list(country_gtype_matrix), fill_color = c("#6a92ff", "#90e25c"), stroke_size = 0, set_name_size = 0, auto_scale = TRUE, text_size=6, show_outside = "none", padding=0.15) + coord_flip()
venn.plot = ggarrange(venn.all, venn.filt, nrow=2, heights=c(1,1))
ggsave(file="../figures/STs_gtype_venn.pdf", width=4, height=8)

# plot novelty
mags.lvs = combined_data[which(combined_data$Genome_Type == "MAG" & !is.na(combined_data$ST_2)),]
mags.lvs$Count = 1
mags.lvs.agg = aggregate(Count ~ ST_2 + Country, data=mags.lvs, FUN=sum)
novel = ggplot(mags.lvs.agg, aes(x=ST_2, y=Count, fill=Country)) +
  geom_bar(stat="identity") +
  theme_classic() +
  scale_fill_manual(values=colorRampPalette(rev(brewer.pal(12, "Set3")))(length(unique(mags.lvs.agg$Country)))) +
  ylab("Number of MAGs") +
  xlab("Number of Locus Variants (LV)") +
  theme(axis.title = element_text(size=18)) +
  theme(axis.text = element_text(size=14)) +
  theme(legend.title = element_text(size=18)) +
  theme(legend.text = element_text(size=14))
ggsave(file="../figures/STs_novel.pdf", width=6, height=6)
