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
  roc.data = map_dfr(.x=outcome, .f=get_sens_spec_lookup, prob_obs)
  roc.data$Genome_type = genome_type
  roc.data$Model = model_name
  return(roc.data)
}

# load best models
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_SamriddhiGupta_Thesis/data/machine_learning/")

# generate roc data
input.files = list.files(path = ".", pattern = ".Rds", recursive = TRUE)
input.files = input.files[which(grepl("inf_filtmags_all|inf_isolates_all", input.files))]
roc.list = lapply(input.files, function(x) {
  cat("Generating ROC data for", x, "...\n")
  genome_type = gsub("_.*", "", basename(x))
  model = strsplit(basename(x), "_")[[1]][[2]]
  roc.data = gen_roc_data(x, "Diseased", genome_type, model)
  roc.data$FPR = 1-roc.data$specificity
  return(roc.data)
})

roc.combined = as.data.frame(rbindlist(roc.list))
roc.combined$FPR = round(roc.combined$FPR, digits=3)
roc.combined$sensitivity = round(roc.combined$sensitivity, digits=3)
roc.agg.mean = aggregate(sensitivity ~ FPR + Genome_type + Model, data=roc.combined, FUN=mean)

# ensure monotonicity
for (gtype in unique(roc.agg.mean$Genome_type)) {
  for (model in unique(roc.agg.mean$Model)) {
    select.rows = which(roc.agg.mean$Genome_type == gtype & roc.agg.mean$Model == model)
    for (n in select.rows[2:length(select.rows)]){
      roc.agg.mean[n,"sensitivity"] = ifelse(roc.agg.mean[n-1,"sensitivity"] > roc.agg.mean[n,"sensitivity"], roc.agg.mean[n-1,"sensitivity"], roc.agg.mean[n,"sensitivity"])
    }
  }
}

# label
roc.agg.mean$Model = recode(roc.agg.mean$Model, "glmnet" = "Ridge Regression", "rf" = "Random Forest", "xgbTree" = "Gradient Boosting")
roc.agg.mean$Label = paste0(roc.agg.mean$Genome_type, ", ", roc.agg.mean$Model)
roc.agg.mean$Label = gsub("mags", "MAGs + Isolates", roc.agg.mean$Label)
roc.agg.mean$Label = gsub("isolates", "Isolates only", roc.agg.mean$Label)


# get AUCs
perf.files = list.files(path = ".", pattern = "performance_results.csv", recursive=TRUE)
perf.files = perf.files[which(grepl("inf_filtmags_all|inf_isolates_all", perf.files))]
for (perf in perf.files) {
  roc.auc = read.csv(perf)
  genome_type = ifelse(grepl("isolates", basename(perf)), "isolates", "mags")
  genome_ren = gsub("isolates", "Isolates only", genome_type)
  genome_ren = gsub("mags", "MAGs + Isolates", genome_ren)
  for (model in unique(roc.auc$method)) {
    auc.value = round(median(roc.auc[which(roc.auc$method == model),"AUC"], na.rm=TRUE),3)
    model.name = recode(model, "glmnet" = "Ridge Regression", "rf" = "Random Forest", "xgbTree" = "Gradient Boosting")
    roc.agg.mean[which(roc.agg.mean$Genome_type == genome_type & roc.agg.mean$Model == model.name),"AUC"] = paste0(genome_ren, ", ", model.name, ", AUROC = ", auc.value)
  }
}

# plot roc curve
roc.curve = ggplot(roc.agg.mean, aes(x=FPR, y=sensitivity, colour=AUC)) +
  geom_line(linewidth=0.5) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
  theme_classic() + 
  theme(panel.grid.minor = element_blank()) + 
  ylim(0,1) +
  xlim(0,1) +
  ylab("True Positive Rate") + 
  xlab("False Positive Rate") +
  scale_colour_manual(values=rev(c("darkolivegreen", "darkolivegreen3", "darkolivegreen2",
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
