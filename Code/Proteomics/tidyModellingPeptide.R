rm(list = ls())
library(tidyverse)
library(LFQService)
library(multcomp)
VIS_PROT <- FALSE

path <- "second_Draft_V2_ContRemoved"
results <- readRDS(file.path(path, "allPreprocessingDataData.Rds"))

# prots <- sample( unique(results$pepIntensityNormalized$protein_Id), 30 )
# pepIntensity <- results$pepIntensityNormalized %>% filter(protein_Id %in% prots)
pepIntensity <- results$pepIntensityNormalized


results$path <- file.path(results$path, "modelling_peptide")
if (!dir.exists(results$path))
{
  dir.create(results$path)
}

results$config_pepIntensityNormalized$table$factorLevel <- 2
pepConfig <- results$config_pepIntensityNormalized


pepConfig$table$factorKeys()
pepConfig$table$hierarchyLevel <- 2
modelName_Interaction  <- "f_Cells_Treatment_Interaction"

formula_interaction <- make_custom_model_lmer("transformedIntensity  ~ class_therapy + Celltype + class_therapy * Celltype + (1 | peptide_Id)",
                                              model_name = modelName_Interaction)

# models_interaction <- model_analyse(pepIntensity,
#                                     formula_interaction,
#                                     modelName_Interaction,
#                                     subject_Id = pepConfig$table$hkeysLevel())
# xx <-model_analyse_summarize(models_interaction$modelProtein,
#                              modelName_Interaction,
#                              subject_Id = pepConfig$table$hkeysLevel())
# visualization <- model_analyse_summarize_vis(xx,pepConfig$table$hkeysLevel())
models_interaction <- workflow_model_analyse(pepIntensity,
                                             formula_interaction,
                                             modelName_Interaction)

models_interaction2 <- models_interaction(results$path)
names(models_interaction2)
mean(models_interaction2$modellingResult$modelProtein$exists_lmer)

mi <- get_complete_model_fit(models_interaction2$modellingResult$modelProtein)

linfct <- linfct_from_model(mi$linear_model[[1]])
lfc <- linfct$linfct_interactions


A <- rbind(
  "p.NO:HSC - c.NO:HSC" = lfc["class_therapyp.NO:CelltypeHSC", ] - lfc["class_therapyc.NO:CelltypeHSC", ],
  "p.HU:HSC - p.NO:HSC" = lfc["class_therapyp.HU:CelltypeHSC", ] - lfc["class_therapyp.NO:CelltypeHSC", ],
  "p.NO:HSC - p.NO:CMP/MEP" = lfc["class_therapyp.NO:CelltypeHSC", ] - lfc["class_therapyp.NO:CelltypeCMP/MEP", ],
  "p.HU:HSC - p.HU:CMP/MEP" = lfc["class_therapyp.HU:CelltypeHSC", ] - lfc["class_therapyp.HU:CelltypeCMP/MEP", ],
  "c.NO:HSC -  c.NO:CMP/MEP" = lfc["class_therapyc.NO:CelltypeHSC", ] -  lfc["class_therapyc.NO:CelltypeCMP/MEP", ],
  "p.NO:CMP/MEP - c.NO:CMP/MEP" = lfc["class_therapyp.NO:CelltypeCMP/MEP", ] - lfc["class_therapyc.NO:CelltypeCMP/MEP", ],
  "p.HU:CMP/MEP - p.NO:CMP/MEP" = lfc["class_therapyp.HU:CelltypeCMP/MEP", ] - lfc["class_therapyp.NO:CelltypeCMP/MEP", ]
)


B <- rbind(
  #"c.NO:HSC - c.NO:CMP/MEP"=lfc["class_therapyc.NO:CelltypeHSC",] - lfc["class_therapyc.NO:CelltypeCMP/MEP",],
  "p.NO:HSC - c.NO:CMP/MEP" = lfc["class_therapyp.NO:CelltypeHSC", ] - lfc["class_therapyc.NO:CelltypeCMP/MEP", ],
  "p.HU:HSC - c.NO:CMP/MEP" = lfc["class_therapyp.HU:CelltypeHSC", ] - lfc["class_therapyc.NO:CelltypeCMP/MEP", ],
  #"p.NO:CMP/MEP - c.NO:CMP/MEP"=lfc["class_therapyp.NO:CelltypeCMP/MEP",] - lfc["class_therapyc.NO:CelltypeCMP/MEP",],
  "p.HU:CMP/MEP - c.NO:CMP/MEP" = lfc["class_therapyp.HU:CelltypeCMP/MEP", ] - lfc["class_therapyc.NO:CelltypeCMP/MEP", ]
)


#CMP/MEP - HSC
#p.NO - c.NO
#(p.NO + p.Prog) - (c.NO + p.HU)

C <- rbind(
  "[p.NO:HSC - p.NO:CMP/MEP] - [c.NO:HSC - c.NO:CMP/MEP]" =
    (lfc["class_therapyp.NO:CelltypeHSC", ] - lfc["class_therapyp.NO:CelltypeCMP/MEP", ]) -
    (lfc["class_therapyc.NO:CelltypeHSC", ] - lfc["class_therapyc.NO:CelltypeCMP/MEP", ]),
  "[p.HU:HSC - p.HU:CMP/MEP] - [p.NO:HSC - p.NO:CMP/MEP]" =
    (lfc["class_therapyp.HU:CelltypeHSC", ] - lfc["class_therapyp.HU:CelltypeCMP/MEP", ]) -
    (lfc["class_therapyp.NO:CelltypeHSC", ] - lfc["class_therapyp.NO:CelltypeCMP/MEP", ]),
  "[p.HU:HSC - p.HU:CMP/MEP] - [c.NO:HSC - c.NO:CMP/MEP]" =
    (lfc["class_therapyp.HU:CelltypeHSC", ] - lfc["class_therapyp.HU:CelltypeCMP/MEP", ]) -
    (lfc["class_therapyc.NO:CelltypeHSC", ] - lfc["class_therapyc.NO:CelltypeCMP/MEP", ])
)

XX <- rbind(A, B, C)
modelProteinF <- models_interaction2$modellingResult$modelProtein %>% dplyr::filter(exists_lmer == TRUE)
contrasts <- workflow_contrasts_linfct(modelProteinF,
                                       modelName = models_interaction2$modellingResult$modelName,
                                       XX,
                                       contrastfun = LFQService::my_contest
)

res <- contrasts(results$path)
