# load libraries
library(ggsankey)
library(ggplot2)
library(dplyr)
library(gplots)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data")
metadata = read.delim("metadata/Metadata_15022025.tsv")
metadata$Health_Status[is.na(metadata$Health_Status)] = "Unknown status"
metadata$Country[is.na(metadata$Country)] = "Unknown country"

# parse names to remove spaces
metadata$Country = gsub("UnitedRepublicofTanzania", "Tanzania", metadata$Country)
metadata$Country = gsub("UnitedStates", "USA", metadata$Country)
metadata$Country = gsub("UnitedKingdom", "UK", metadata$Country)
metadata$Country = gsub("ElSalvador", "El Salvador", metadata$Country)
metadata$Country = gsub("SouthAfrica", "South Africa", metadata$Country)

# parse data
df = metadata %>%
  make_long(Genome_Type, Health_Status, Country)
df

# tally groups
dagg = df %>%
  dplyr::group_by(node)%>%
  tally()
df2 = merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)


# define colors
sankey.cols = c("#90e25c", "#6a92ff", "#0072B2", "#D55E00", "grey70")
names(sankey.cols) = c("MAG", "Isolate", "Healthy", "Diseased", "Unknown status")

# plot
ggplot(df2, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
               label = paste0(node," n=", n))) +
  geom_sankey(flow.alpha = 0.5, node.color = 0, alpha=0.6, type="sankey", space=18) +
  scale_fill_manual(values=sankey.cols) +
  geom_sankey_label(size = 4, color = 1, fill = "white", alpha=0.75, label.size = NA, space=18) +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none") +
  theme(axis.title = element_blank()) +
  theme(axis.text = element_blank())
ggsave(file="figures/metadata_sankey.pdf", height=9, width=11)
