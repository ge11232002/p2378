rm(list = ls())
library(tidyverse)
library(LFQService)
library(multcomp)
VIS_PROT <- FALSE

outpath <- "second_Draft_V4_WithContaminants"
proteins_fun <- readRDS(file = file.path(outpath, "protintensity_fun.Rda"))

outpath <- "modelling_results"

config <- proteins_fun("unnest")$config
config$table$factorKeys()
config$table$fkeysLevel()
config$table$workIntensity
#config$table$factorLevel <- 2 
#config$table$hierarchyLevel <- 1


protIntensity <- proteins_fun("unnest")
protIntensityNormalized <- protIntensity$data
config <- protIntensity$config

modelling_dir <- "protein_modelling_results_v1"


config$table$fkeysLevel()
as.character(unique(LFQService::make_interaction_column_config(protIntensityNormalized, config,sep = ":")$interaction))


modelName  <- "f_Cells_Treatment_Interaction_Protein"
model <- "~ celltype_ + class_therapy_progressed_ + class_therapy_progressed_ * celltype_"


model <- paste0(config$table$getWorkIntensity() ,model)
modelFunction <- make_custom_model_lm( model, model_name = "Model")


reportColumns = c("moderated.p.value",
                  "moderated.p.value.adjusted")

DEBUG <- FALSE

Contrasts <- c(
  "(patient.NO.yes + patient.HU.yes)/2" = "(class_therapy_progressed_patient.NO.yes + class_therapy_progressed_patient.HU.yes)/2",
  "(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2)/2" = "(class_therapy_progressed_patient.NO.no + `(patient.NO.yes + patient.HU.yes)/2`)/2",
  "(control.NO.no + patient.HU.no)/2" = "(class_therapy_progressed_control.NO.no + class_therapy_progressed_patient.HU.no)/2",
  
  "PV_progressed and PV_NO vs Control and PV_HU" = "`(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2)/2` - `(control.NO.no + patient.HU.no)/2`",
  
  "celltype_CMP/MEP:(patient.NO.yes + patient.HU.yes)/2" = "(`celltype_CMP/MEP:class_therapy_progressed_patient.NO.yes` + `celltype_CMP/MEP:class_therapy_progressed_patient.HU.yes`)/2",
  "celltype_CMP/MEP:(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2)/2" = "(`celltype_CMP/MEP:class_therapy_progressed_patient.NO.no` + `celltype_CMP/MEP:(patient.NO.yes + patient.HU.yes)/2`)/2",
  "celltype_CMP/MEP:(control.NO.no + patient.HU.no)/2" = "(`celltype_CMP/MEP:class_therapy_progressed_control.NO.no` + `celltype_CMP/MEP:class_therapy_progressed_patient.HU.no`)/2",
  
  "CMP/MEP:PV_progressed and PV_NO vs Control and PV_HU" = "`celltype_CMP/MEP:(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2)/2` - `celltype_CMP/MEP:(control.NO.no + patient.HU.no)/2`",
  
  "celltype_HSC:(patient.NO.yes + patient.HU.yes)/2" = "(`celltype_HSC:class_therapy_progressed_patient.NO.yes` + `celltype_HSC:class_therapy_progressed_patient.HU.yes`)/2",
  "celltype_HSC:(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2)/2" = "(`celltype_HSC:class_therapy_progressed_patient.NO.no` + `celltype_HSC:(patient.NO.yes + patient.HU.yes)/2`)/2",
  "celltype_HSC:(control.NO.no + patient.HU.no)/2" = "(`celltype_HSC:class_therapy_progressed_control.NO.no` + `celltype_HSC:class_therapy_progressed_patient.HU.no`)/2",
  
  "HSC:PV_progressed and PV_NO vs Control and PV_HU" = "`celltype_HSC:(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2)/2` - `celltype_HSC:(control.NO.no + patient.HU.no)/2`",
  
  "differences : HSC vs CMP/MEP for progression" = "`HSC:PV_progressed and PV_NO vs Control and PV_HU` - `CMP/MEP:PV_progressed and PV_NO vs Control and PV_HU`"
)



res <- application_run_modelling_V2(outpath = outpath,
                                    data = protIntensityNormalized,
                                    pepConfig = config,
                                    modelFunction = modelFunction,
                                    Contrasts = Contrasts,
                                    modelling_dir =modelling_dir)
