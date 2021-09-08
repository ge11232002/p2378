rm(list=ls())
library(tidyverse)
library(LFQService)
library(multcomp)
VIS_PROT <- FALSE

path <- "second_Draft_V2"
results <- readRDS(file.path(path,"allPreprocessingDataData.Rds"))
prots <- sample( unique(results$pepIntensityNormalized$protein_Id), 30 )
pepIntensity <- results$pepIntensityNormalized %>% filter(protein_Id %in% prots)

dim(results$pepIntensityNormalized)
dim(pepIntensity)
results$path <- file.path(results$path, "modelling")
if(!dir.exists(results$path))
{
  dir.create(results$path)
}

results$config_pepIntensityNormalized$table$factorLevel <- 2
pepConfig<- results$config_pepIntensityNormalized
pepConfig$table$factorKeys()

modelName_Interaction  <- "f_Cells_Treatment_Interaction"
formula_interaction <- make_custom_model_lmer("transformedIntensity  ~ class_therapy + Celltype + class_therapy * Celltype + (1 | peptide_Id)")
models_interaction <- workflow_model_analyse(pepIntensity, formula_interaction, modelName_Interaction)
models_interaction <- models_interaction(results$path)


mi <- get_complete_model_fit(models_interaction$summaryResult$modelProteinF)

linfct <- linfct_from_model(mi)
rownames(linfct$linfct_interactions)

mk_group_differences <- function(linfct_factors){
  group_differences <- rbind(
    "Celltype: HSC - CMP/MEP" = linfct_factors["CelltypeHSC",] - linfct_factors["CelltypeCMP/MEP",],
    "class_therapy: p.NO - c.NO" = linfct_factors["class_therapyp.NO",] - linfct_factors["class_therapyc.NO",],
    "class_therapy: p.HU - c.NO" = linfct_factors["class_therapyp.HU",] - linfct_factors["class_therapyc.NO",],
    "class_therapy: p.HU - p.NO" = linfct_factors["class_therapyp.HU",] - linfct_factors["class_therapyp.NO",]
  )
  return(group_differences)
}

XX <- rbind(linfct_factors_contrasts(mi) , linfct_all_possible_contrasts(linfct$linfct_interactions))
contrasts <- workflow_contrasts_linfct_ALL(models_interaction$summaryResult,XX)
contrasts(results$path)
