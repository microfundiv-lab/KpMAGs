# load libraries
library(mikropml)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

# function to get tpr and fpr per outcome
get_sens_spec_lookup = function(outcome, data){
  selected = data %>%
    select(!!as.name(outcome), observed)
  
  total = selected %>%
    count(observed) %>%
    pivot_wider(names_from=observed, values_from=n)
  
  selected %>%
    arrange(desc(!!as.name(outcome))) %>%
    mutate(is_outcome = (observed == outcome),
           tp = cumsum(is_outcome),
           fp = cumsum(!is_outcome),
           precision = tp / (tp + fp),
           sensitivity = tp / as.numeric(total[outcome]),
           fpr = fp / (sum(total)-as.numeric(total[outcome])),
           specificity = 1-fpr) %>% # TN/TN+FP = 1-FPR
    add_column(class = outcome) %>%
    select(precision, sensitivity, specificity, fpr, class)
}

# function to generate roc data per model
gen_roc_data = function(model, outcome, genome_type, model_name){
  select.model = readRDS(model)
  prob = predict(select.model$trained_model, select.model$test_data, type="prob")
  observed = select.model$test_data$Variable
  prob_obs = bind_cols(prob, observed=observed)
  roc.data = purrr::map_dfr(.x=outcome, .f=get_sens_spec_lookup, prob_obs)
  roc.data$Genome_type = genome_type
  roc.data$Model = model_name
  return(roc.data)
}

# load best models
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/machine_learning/models")

# generate roc data
input.files = list.files(path = ".", pattern = ".Rds", recursive = FALSE)
roc.list = lapply(input.files, function(x) {
  cat("Generating ROC data for", x, "...\n")
  genome_type = gsub("_.*", "", basename(x))
  model = strsplit(basename(x), "_")[[1]][[2]]
  roc.data = gen_roc_data(x, "Diseased", genome_type, model)
  roc.data$FPR = 1-roc.data$specificity
  return(roc.data)
})

roc.combined = as.data.frame(rbindlist(roc.list))
roc.combined$FPR = round(roc.combined$FPR, digits=4)
roc.combined$sensitivity = round(roc.combined$sensitivity, digits=4)
roc.agg.median = aggregate(sensitivity ~ FPR + Genome_type + Model, data=roc.combined, FUN=median)

# ensure monotonicity
for (gtype in unique(roc.agg.median$Genome_type)) {
  for (model in unique(roc.agg.median$Model)) {
    select.rows = which(roc.agg.median$Genome_type == gtype & roc.agg.median$Model == model)
    for (n in select.rows[2:length(select.rows)]){
      roc.agg.median[n,"sensitivity"] = ifelse(roc.agg.median[n-1,"sensitivity"] > roc.agg.median[n,"sensitivity"], roc.agg.median[n-1,"sensitivity"], roc.agg.median[n,"sensitivity"])
    }
  }
}

# label
roc.agg.median$Model = recode(roc.agg.median$Model, "glmnet" = "Ridge Regression", "rf" = "Random Forest", "xgbTree" = "Gradient Boosting")
roc.agg.median$Label = paste0(roc.agg.median$Genome_type, ", ", roc.agg.median$Model)
roc.agg.median$Label = gsub("mags", "MAGs + Isolates", roc.agg.median$Label)
roc.agg.median$Label = gsub("isolates", "Isolates only", roc.agg.median$Label)


# get AUCs
perf.files = list.files(path = "../", pattern = "_results.csv", recursive=FALSE, full.names = TRUE)
for (perf in perf.files) {
  roc.auc = read.csv(perf)
  genome_type = ifelse(grepl("isolates", basename(perf)), "isolates", "mags")
  genome_ren = gsub("isolates", "Isolates only", genome_type)
  genome_ren = gsub("mags", "MAGs + Isolates", genome_ren)
  for (model in unique(roc.auc$method)) {
    auc.value = round(median(roc.auc[which(roc.auc$method == model),"AUC"], na.rm=TRUE),3)
    model.name = recode(model, "glmnet" = "Ridge Regression", "rf" = "Random Forest", "xgbTree" = "Gradient Boosting")
    roc.agg.median[which(roc.agg.median$Genome_type == genome_type & roc.agg.median$Model == model.name),"AUC"] = paste0(genome_ren, ", ", model.name, ", AUROC = ", auc.value)
  }
}

# plot roc curve
roc.curve = ggplot(roc.agg.median, aes(x=FPR, y=sensitivity, colour=AUC)) +
  geom_line(linewidth=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
  theme_classic() + 
  theme(panel.grid.minor = element_blank()) + 
  ylim(0,1) +
  xlim(0,1) +
  ylab("True Positive Rate") + 
  xlab("False Positive Rate") +
  scale_colour_manual(values=rev(c("darkolivegreen2", "darkolivegreen3", "darkolivegreen",
                               "tomato1", "tomato3", "tomato4"))) +
  theme(legend.position="right") +
  guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(strip.background = element_blank(), strip.text = element_text(size=14)) +
  theme(legend.position=c(0.6,0.2), legend.box = "horizontal", legend.text=element_text(size=12),
        legend.title = element_blank()) +
  theme(axis.title.y = element_text(size=16)) + 
  theme(axis.text.y = element_text(size=14)) + 
  theme(axis.title.x = element_text(size=16)) + 
  theme(axis.text.x = element_text(size=14))
