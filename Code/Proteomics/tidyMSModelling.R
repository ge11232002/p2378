rm(list=ls())
library(tidyverse)
library(LFQService)
library(multcomp)
VIS_PROT <- FALSE

results <- readRDS("allData.Rds")

results$pepIntensityNormalized
results$config_pepIntensityNormalized$table$factorLevel <- 2
pepConfig<- results$config_pepIntensityNormalized
pepConfig$table$factorKeys()

pepConfig$table$factorKeys()

modelName  <- "f_ClassTherapy_Celltype_r_Patient"
if(TRUE){
  formula_randomPeptide <- make_custom_model_lmer("transformedIntensity  ~ class_therapy + Celltype  + (1|Patient_ID) + (1 | peptide_Id)")
  res_model_REPEATED <- workflow_no_interaction_modelling(results$pepIntensityNormalized,
                                                          results$config_pepIntensityNormalized,
                                                          formula_randomPeptide, modelName)
  saveRDS(res_model_REPEATED,file=file.path(results$path, paste0(modelName,".rda")))
}
res_model_REPEATED <- readRDS(file.path(results$path, paste0(modelName,".rda"))) 


modelName  <- "f_ClassTherapy_Celltype"
if(TRUE){
  formula_randomPatient <- make_custom_model_lmer("transformedIntensity  ~ class_therapy + Celltype  + (1|peptide_Id)")
  res_model <- workflow_no_interaction_modelling(results,formula_randomPatient, modelName)
  saveRDS(res_model,file=file.path(results$path, paste0(modelName,".rda")))
}

res_model <- readRDS(file.path(results$path, paste0(modelName,".rda"))) 


lltest <- workflow_likelihood_ratio_test(res_model$modelProteinF ,
                                         res_model$modelName ,
                                         res_model_REPEATED$modelProteinF ,
                                         res_model_REPEATED$modelName )

pdf(file.path(results$path, "LRT_histogram.pdf"))
hist(lltest$likelihood_ratio_test.pValue, breaks=20)
dev.off()

res_model$TwoFactorModelFactor2 <- inner_join(res_model$TwoFactorModelFactor2,lltest , by = "protein_Id")
res_model_REPEATED$TwoFactorModelFactor2 <- inner_join(res_model_REPEATED$TwoFactorModelFactor2,lltest , by = "protein_Id")

readr::write_csv(res_model$TwoFactorModelFactor2, path = file.path(results$path, paste0("Contrasts_Auto_SignificanceValues_",  res_model$modelName, ".csv")))
readr::write_csv(res_model_REPEATED$TwoFactorModelFactor2, path = file.path(results$path, paste0("Contrasts_Auto_SignificanceValues_",  res_model_REPEATED$modelName, ".csv")))

res_cond_r_pep_Pivot <- pivot_model_contrasts_2_Wide(res_model$TwoFactorModelFactor2)
write_csv(res_cond_r_pep_Pivot, path=file.path(results$path, paste0("Contrasts_SignificanceValues_", res_model$modelName, "_PIVOT.csv")))

res_cond_r_pep_r_pat_Pivot <- pivot_model_contrasts_2_Wide(res_model_REPEATED$TwoFactorModelFactor2)
write_csv(res_cond_r_pep_r_pat_Pivot, path=file.path(results$path, paste0("Contrasts_SignificanceValues_", res_model_REPEATED$modelName, "_PIVOT.csv")))

tmp <- inner_join(res_model$TwoFactorModelFactor2,  res_model_REPEATED$TwoFactorModelFactor2, by=c("protein_Id", "factor", 'isSingular', 'peptide_Id_n', 'lhs'))

pdf(file.path(results$path, "compareModels.pdf"))
ggplot(tmp, aes(x = estimate.x, y = estimate.y, )) + geom_point() + facet_wrap(~lhs) 
ggplot(tmp, aes(x = p.value.x, y = p.value.y, )) + geom_point() + facet_wrap(~lhs) 
dev.off()


# add group averages -----

m <- res_model$modelProteinF$lmer_f_ClassTherapy_Celltype[[1]]
linfct <- lmer4_linfct_from_model(m)

res_cond_r_pep_grA <- workflow_group_averages(res_model$modelProteinF,res_model$modelName,  results$path, linfct$linfct_interactions, lltest )
res_cond_r_pep_r_pat_grA <- workflow_group_averages(res_model_REPEATED$modelProteinF,res_model_REPEATED$modelName,  results$path, linfct$linfct_interactions,lltest )

write_csv(res_cond_r_pep_grA, path=file.path(results$path, paste0("GroupAverages_", res_model$modelName, ".csv")))
write_csv(res_cond_r_pep_r_pat_grA, path=file.path(results$path, paste0("GroupAverages_", res_model_REPEATED$modelName, ".csv")))


# Doing protein plots with prediction -----
if(FALSE){
  
  library(lme4)
  res_model$modelProteinF <- res_model$modelProteinF %>%
    mutate(nRandom_Interaction_plot = map2(!!sym(paste0("lmer_", res_model$modelName)), protein_Id, plot_model_and_data))
  
  
  #plot_model_and_data_TWO(res_model_REPEATED$modelProteinF$lmer_f_Condition_r_peptid_r_patient[[1]],"A",legend.position = "bottom", firstlast = TRUE)
  #plot_model_and_data_TWO(res_model_REPEATED$modelProteinF$lmer_f_Condition_r_peptid_r_patient[[1]],"A",legend.position = "bottom", firstlast = FALSE)
  
  res_model_REPEATED$modelProteinF <- res_model_REPEATED$modelProteinF %>%
    mutate(nRandom_Interaction_plot = map2(!!sym(paste0("lmer_", res_model_REPEATED$modelName)), protein_Id, plot_model_and_data_TWO, firstlast = TRUE))
  
  tmp <- inner_join(dplyr::select(res_model$modelProteinF, protein_Id, nRandom_Interaction_plot),
                    dplyr::select(res_model_REPEATED$modelProteinF, protein_Id, nRandom_Interaction_plot), by="protein_Id" )
  
  tmp <- inner_join(tmp, lltest)
  
  tmp <- tmp %>% filter(likelihood_ratio_test.pValue < 0.000001)
  pdf(file.path(results$path,"Significant_Random_Patient.pdf"), width =8, height =8)
  for(i in 1:nrow(tmp)){
    print(gridExtra::grid.arrange(tmp$nRandom_Interaction_plot.x[[i]], tmp$nRandom_Interaction_plot.y[[i]]))
  }
  dev.off()
  
}