# load libraries
library(ggpubr)

# run external scripts
source("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/scripts/alex/ml_mags-vs-isolates.R")
source("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/scripts/alex/pyseer_volcano-cog.R")

# combine figure
ggarrange(gwas.fig, ml.fig, ncol=1)
ggsave("../figures/gwas-ml_combined.pdf", height=10, width=14)
