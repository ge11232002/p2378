####################### End of configuration #################################


# create result structure
modelling_path <- file.path(outpath, modelling_dir)


if(!dir.exists(outpath)){
  dir.create(outpath)
}


if(!dir.exists(modelling_path)){
  dir.create(modelling_path)
}


###  Just run the modelling


config$table$hierarchyLevel <- 2
wideFRAME <- LFQService::toWideConfig(protIntensityNormalized,
                                      config)

lfq_write_table(separate_hierarchy(protIntensityNormalized,
                                   config),
                file.path(modelling_path, "protein_intensities.csv"))


# render protein quantification reprot
if(DEBUG){
LFQService::render_MQSummary_rmd(protIntensityNormalized,
                                 config$clone(deep=TRUE),
                                 pep=FALSE,
                                 workdir = ".",
                                 dest_path = modelling_path,
                                 dest_file_name="protein_intensities_qc.pdf")
}

#################################################
### Do missing value imputation

res_contrasts_imputed <- workflow_missigness_impute_contrasts(protIntensityNormalized,
                                                              config,
                                                              Contrasts)
## linear models ----
pepConfig <- config

# first model ----

model <- paste0(pepConfig$table$getWorkIntensity() , model)

modelFunction <- make_custom_model_lm( model , modelName )
modelFunction$model_fun(get_formula = T)

### make contrasts -----

modellingResult_fun <- workflow_model_analyse(protIntensityNormalized,
                                              modelFunction,
                                              modelName,
                                              subject_Id = pepConfig$table$hkeysLevel())

modellingResult <- modellingResult_fun()

#names(modellingResult)

m <- get_complete_model_fit(modellingResult$modellingResult$modelProtein)

#factor_contrasts <- linfct_factors_contrasts(m$linear_model[[1]])
linfct <- linfct_from_model(m$linear_model[[1]], as_list = FALSE)
linfct_A <- linfct_matrix_contrasts(linfct, Contrasts)

if(exists("relevant_contrasts")){ 
  linfct_A <- linfct_A[relevant_contrasts ,, drop=F]
}

if(DEBUG){
  pdf(file.path(modelling_path,"Linear_functions.pdf"), width=18, height=10)
  quantable::imageWithLabels(t(linfct_A),
                             col = quantable::getBlueWhiteRed(),
                             marLeft = c(8,10,4.1,2.1))
  dev.off()
}

#rownames(factor_contrasts$linfct_interactions)

modelProteinF <- modellingResult$modellingResult$modelProtein %>% dplyr::filter(exists_lmer == TRUE)

res_contrasts <- workflow_contrasts_linfct(modelProteinF ,
                                           modellingResult$modellingResult$modelName,
                                           linfct_A,
                                           pepConfig,
                                           prefix =  "Contrasts",
                                           contrastfun = modelFunction$contrast_fun)

#names(xx)
#View(xx$contrast_result)

xx <- res_contrasts(modelling_path,columns = reportColumns)
xx_imputed <- res_contrasts_imputed("long",what = "contrasts")

merge_contrasts_results <- function(contrast_minimal,
                                    xx_imputed, subject_Id){
  res <- right_join(contrast_minimal, xx_imputed, by=c(subject_Id,"lhs" = "contrast"))
  res <- res %>% rename(contrast = lhs)
  res <- res %>% mutate(pseudo_estimate = case_when(is.na(estimate) ~ imputed, TRUE ~ estimate))
  res <- res %>% mutate(is_pseudo_estimate = case_when(is.na(estimate) ~ TRUE, TRUE ~ FALSE))
  
  for(column in reportColumns){
    print(column)
    res <- res %>% mutate(!!sym(paste0("pseudo_",column)) := case_when(is.na(estimate) ~ 0, TRUE ~ !!sym(column)))
  }
  res <- res %>% dplyr::select(-imputed, -meanArea)
  return(res)
}

contrast_results <- merge_contrasts_results(xx$contrast_minimal, xx_imputed,
                                            subject_Id = pepConfig$table$hkeysLevel())

separate_hierarchy(contrast_results, config) -> contrast_results
lfq_write_table(contrast_results, path = file.path(modelling_path, "foldchange_estimates.csv"))
