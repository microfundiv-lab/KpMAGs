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

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/metadata")
combined_data = read.delim("Metadata_09122024.tsv")
combined_data = combined_data[which(combined_data$ST != "0"),]

# subset countries
country_health_table = table(combined_data$Country, combined_data$Health_Status)
country_health_matrix = as.data.frame.matrix(country_health_table)
select.countries = rownames(country_health_matrix)[which(country_health_matrix$Diseased > 1 & country_health_matrix$Healthy > 1 )]
combined_data_filt = combined_data[which(combined_data$Country %in% select.countries),]
country_health_table = table(combined_data_filt$ST, combined_data_filt$Health_Status)
country_health_matrix = as.data.frame.matrix(country_health_table)

# ST analysis
ST_health_table = table(combined_data$ST, combined_data$Health_Status)
ST_health_matrix = as.data.frame.matrix(ST_health_table)
ST_health_df = as.data.frame(ST_health_table)

# check diversity
dis.div = diversity(ST_health_matrix$Diseased)
health.div = diversity(ST_health_matrix$Healthy)

# extract top
top_diseased_STs = combined_data %>%
  filter(Health_Status == "Diseased") %>%
  group_by(ST) %>%
  summarize(Count = n()) %>%
  arrange(desc(Count)) %>%
  top_n(6, Count)

top_healthy_STs = combined_data %>%
  filter(Health_Status == "Healthy") %>%
  group_by(ST) %>%
  summarize(Count = n()) %>%
  arrange(desc(Count)) %>%
  top_n(5, Count)

#plot for diseased Population
plot_diseased = ggplot(top_diseased_STs, aes(x = reorder(ST, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = '#D55E00', width = 0.4) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 22, margin = margin(r = 10)),  
    axis.title.x = element_text(size = 22), 
    axis.text.x = element_text(size = 20), 
    axis.text.y = element_text(size = 20),
    plot.margin=unit(c(1.5,1,1,1), 'cm'),
    plot.title = element_blank(),  
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dashed'),
    panel.grid.minor = element_blank() 
  ) 

#healthy population
plot_healthy = ggplot(top_healthy_STs, aes(x = reorder(ST, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = '#0072B2', width = 0.4) +
  labs(x = "Sequence Type (ST)", y = "") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 22, margin = margin(r = 10)),  
    axis.title.x = element_text(size = 22),  
    axis.text.x = element_text(size = 20),  
    axis.text.y = element_text(size = 20),  
    plot.title = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, linetype = 'dashed'),
    panel.grid.minor = element_blank() 
  )

combined_plot = ggarrange(plot_diseased, plot_healthy, align="v", ncol=1, heights=c(1.25,1))

# get values for venn diagram
gen_venn_list = function(x) {
  disease.list = rownames(x)[which(x$Diseased > 0)]
  health.list = rownames(x)[which(x$Healthy > 0)]
  venn.list = list(Diseased = disease.list, Healthy = health.list)
  return(venn.list)
}

venn.all = ggvenn(gen_venn_list(ST_health_matrix), fill_color = c("#D55E00", "#0072B2"), stroke_size = 0, set_name_size = 7, auto_scale = TRUE, text_size=6, show_outside = "none", padding=0.15)
venn.filt = ggvenn(gen_venn_list(country_health_matrix), fill_color = c("#D55E00", "#0072B2"), stroke_size = 0, set_name_size = 0, auto_scale = TRUE, text_size=6, show_outside = "none", padding=0.15)
venn.plot = ggarrange(venn.all, venn.filt, nrow=2, heights=c(1.45,1))

# arrange
arrang.plot = ggarrange(combined_plot, venn.plot, ncol=1, font.label = list(size = 24, color = "black", face = "bold", family = NULL), widths=c(1.3,1))
annotate_figure(arrang.plot, left = textGrob("Number of genomes", rot = 90, vjust = 3, hjust=-0.5, gp = gpar(cex = 2)))
ggsave(file="../figures/STs.pdf", width=9, height=14)
