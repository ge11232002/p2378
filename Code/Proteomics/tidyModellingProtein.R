rm(list=ls())
library(tidyverse)
library(LFQService)
library(multcomp)
VIS_PROT <- FALSE

path <- "second_Draft_V4_WithContaminants"
proteins_fun <- readRDS(file = file.path(path, "protintensity_fun.Rda"))



path <- file.path(path, "modelling_protein")
if(!dir.exists(path))
{
  dir.create(path)
}

pepConfig<- proteins_fun("unnest")$config
pepConfig$table$factorKeys()
pepConfig$table$fkeysLevel()
pepConfig$table$workIntensity
protIntensity <- LFQService::applyToIntensityMatrix(proteins_fun("unnest")$data, pepConfig, .func = robust_scale) 

pepConfig$table$fkeysLevel()
as.character(unique(LFQService::make_interaction_column_config(protIntensity, pepConfig,sep = ":")$interaction))



if(FALSE){
patient.NO.no = PV_NO
control.NO.no  = Control
patient.HU.no = PV_HU
patient.NO.yes = PV_progressed
patient.HU.yes = PV_progressed
}

#(PV_NO + PV_progressed) - (Control + PV_HU)

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


#(patient.NO.no + (patient.NO.yes + patient.HU.yes)/2) - (control.NO.no + patient.HU.no)

#patient_NO_no - control_NO_no
#patient_HU_no - patient_NO_no
#patient_HU_no - control_NO_no

#(patient_HU_yes + patient_NO_yes)/2 - patient_NO_no
#(patient_HU_yes + patient_NO_yes)/2 - patient_HU_no
#(patient_HU_yes + patient_NO_yes)/2 - control_NO_no




# linear modelling

modelName_Interaction  <- "f_Cells_Treatment_Interaction_Protein"
formula_interaction <- make_custom_model_lm("medpolish  ~ celltype_ + class_therapy_progressed_ + class_therapy_progressed_ * celltype_", modelName_Interaction)
models_interaction <- workflow_model_analyse(protIntensity, formula_interaction, modelName_Interaction)
#dd <- models_interaction()
models_interaction <- models_interaction()



m <- get_complete_model_fit(models_interaction$modellingResult$modelProtein)

linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
linfct_A <- linfct_matrix_contrasts(linfct, Contrasts)
rownames(linfct_A)

linfct_A <- linfct_A[c(4,8,12,13),, drop=F]

par(mar=c( c(8,15,4.1,2.1)))
quantable::imageWithLabelsNoLayout(t(linfct_A), col = quantable::getBlueScale())
dev.off()

modelProteinF <- models_interaction$modellingResult$modelProtein %>% dplyr::filter(exists_lmer == TRUE)

res <- workflow_contrasts_linfct(modelProteinF ,
                                 models_interaction$modellingResult$modelName,
                                 linfct_A,
                                 pepConfig,
                                 prefix =  "Contrasts",
                                 contrastfun = LFQService::my_contrast_V2)






xx <- res(path,columns = c("moderated.p.value",
                           "moderated.p.value.adjusted"))

xx <- res()


xx$visualization$Contrasts_Histogram_moderated.p.value
xx$visualization$Contrasts_Volcano_moderated.p.value
xx$visualization$Contrasts_Volcano_p.value
xx$visualization$Contrasts_Histogram_p.value

View(res_contrast)
res_contrast <- workflow_missigness_impute_contrasts(protIntensity, pepConfig, Contrasts)
xxx_imputed <- res_contrast("wide",what = "all")
xxx_imputed <- res_contrast("wide",what = "factors")
dim(xxx_imputed)
xxx_imputed_contrasts <- res_contrast("raw",what='contrasts')
xxx_imputed_contrasts <- xxx_imputed_contrasts %>% filter(contrast %in% rownames(linfct_A))

dd <- xxx_imputed_contrasts %>% unite(contrast.v , value, contrast, sep="~") %>% spread(contrast.v, int_val)


r <- left_join(dd,xx$contrasts_wide)
lfq_write_table(separate_hierarchy(r, pepConfig), path=file.path(path, "AllData_WideFormat.csv"))

xxx_imputed_factors <- res_contrast("wide",what='factors')
lfq_write_table(xxx_imputed_factors, path=file.path(path, "Interaction_Intensities.csv"))






