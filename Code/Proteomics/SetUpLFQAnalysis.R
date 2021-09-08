rm(list = ls())

library(readr)
library(tidyverse)
library(LFQService)
library(tidyr)
library(dplyr)

HEATMAP <- TRUE
ALLPROTEINPLOTS <- FALSE
MQSUMMARY <- FALSE

path <- "second_Draft_V4_WithContaminants"
path_qc <- file.path(path, "path_qc")


if(!dir.exists(path)){
  dir.create(path)
}

if(!dir.exists(path_qc)){
  dir.create(path_qc)
}


resPepProtAnnot <- read_rds(path = "data/20180923_230031_SA_1807FMA_HSPCPV5frwoisolw_SNP_Report.rds")
(table(stringr::str_count(resPepProtAnnot$PG.ProteinAccessions, pattern = ";")))

createSpectronautPeptideConfiguration <- function(isotopeLabel="isotope",
                                                  qValue="EG.Qvalue"){
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "filename"
  
  # measurement levels.
  atable$hierarchy[["protein_Id"]] <- "PG.ProteinAccessions"
  atable$hierarchy[["peptide_Id"]] <- "PEP.StrippedSequence"
  atable$hierarchy[["modPeptide_Id"]] <- "EG.ModifiedSequence"
  atable$hierarchy[["precursor_Id"]] <- c("EG.ModifiedSequence","FG.Charge")#EG.PrecursorId
  
  #
  atable$ident_qValue = qValue
  atable$workIntensity = "FG.Quantity"
  #atable$workIntensity = "FG.MS1PeakArea"
  atable$isotopeLabel = isotopeLabel
  
  atable$factors[["celltype_"]] = "celltype"
  atable$factors[["class_therapy_progressed_"]] = "class_therapy_progressed"
  atable$factors[["Patient_ID"]] = "id"
  
  atable$factorDepth = 2
  anaparam <- AnalysisParameters$new()
  anaparam$min_peptides_protein <- 2
  configuration <- AnalysisConfiguration$new(atable, anaparam)
  return(configuration)
}


complete_cases <- function(data, config){
  message("completing cases")
  fkeys <- c(config$table$fileName , config$table$sampleName, config$table$factorKeys())
  hkeys <- c(config$table$isotopeLabel, config$table$hierarchyKeys())
  res <- data %>% tidyr::complete(
    !!sym(config$table$fileName),
    tidyr::nesting(!!!syms(hkeys))
  )
  return(res)
}

config <- createSpectronautPeptideConfiguration()
precursorData <- setup_analysis(resPepProtAnnot, config, cc = FALSE)
dim(complete_cases(precursorData, config))

precursorData <- remove_small_intensities(precursorData, config,threshold = 4)
# filter qvalues and aggregate peptides -----


#config$parameter$maxqVal_experiment_threshold

hierarchy_counts(precursorData,config)
config$parameter$qVal_experiment_threshold = 0.001
config$parameter$qVal_individual_threshold = 0.01

data_NA_QVal <- filter_byQValue(precursorData, config)
hierarchy_counts(data_NA_QVal,config)

#data_NA_QVal %>% dplyr::select(c(config$table$idVars(),config$table$valueVars()))
#setdiff(config$table$idVars(), c("modPeptide_Id","precursor_Id"))

resDataStart <- data_NA_QVal %>% group_by_at(setdiff(config$table$idVars(), c("modPeptide_Id","precursor_Id"))) %>%
  summarize(!!config$table$getWorkIntensity() := sum(!!sym(config$table$getWorkIntensity()), na.rm = TRUE),
            EG.Qvalue = min(EG.Qvalue,na.rm = TRUE)) %>% ungroup()

pepConfig <- config$clone(deep = TRUE)
pepConfig$table$hierarchy$modPeptide_Id <- NULL
pepConfig$table$hierarchy$precursor_Id <- NULL
config <- pepConfig
config$table$hierarchyDepth <- 1
config$table$hkeysDepth()

# This code should be the same for maxquant ----
#resDataStart <- LFQService::make_interaction_column_config(resDataStart, config)

readr::write_csv(resDataStart, path = file.path(path_qc, "annotatedTable_Peptide_RAW_Data.csv"))


if (MQSUMMARY) {
  LFQService::render_MQSummary_rmd(resDataStart, config$clone(deep = TRUE) , dest_path = path_qc,  workdir = ".")
}


config$table$workIntensity <- config$table$workIntensity[1]
config$table$is_intensity_transformed <- FALSE
x3 <- summarize_hierarchy(resDataStart, config)
#x3 %>% inner_join(resDataStart, by="protein_Id") -> resDataStart

############################
# Start filtering
hierarchy_counts(resDataStart,config)

results <- LFQService::workflow_MQ_protoV1(resDataStart,
                               config,
                               path,
                               peptideFilterFunction = LFQService:::.workflow_MQ_filter_peptides_V3 )

hierarchy_counts(results$pepIntensityNormalized,results$config_pepIntensityNormalized)

saveRDS(results, file = file.path(path, "peptideLevel.Rda"))


# results$nrPeptidesPerProtein_filtered
# results$nrPeptidesPerProtein_start
summarize_hierarchy(results$pepIntensityNormalized,results$config_pepIntensityNormalized)

############# Quantify Proteins ################
results$pepIntensityNormalized <- LFQService::complete_cases(results$pepIntensityNormalized, results$config_pepIntensityNormalized)
protintensity <- medpolish_protein_quants( results$pepIntensityNormalized, results$config_pepIntensityNormalized )
saveRDS(protintensity, file = file.path(path, "protintensity_fun.Rda"))


### READ data ----
library(LFQService)
library(tidyverse)
path <- "second_Draft_V4_WithContaminants"
path_qc <- file.path(path, "path_qc")

protintensity <- readRDS(file = file.path(path, "protintensity_fun.Rda"))

lfq_write_table(LFQService::separate_hierarchy(protintensity("unnest")$data, config),
                 path = file.path(path_qc,"transformed_ProteinIntensities.csv"))

lfq_write_table(LFQService::separate_hierarchy(protintensity("wide")$data, config), path = file.path(path_qc,"transformed_ProteinIntensities_WIDE.csv"))


pdf(file.path(path_qc,"heatmapSampleCorProteinLevel.pdf"), width = 10, height = 10)
ggg <- plot_heatmap_cor(protintensity("unnest")$data, protintensity("unnest")$config)
print(ggg)
dev.off()

pdf(file.path(path_qc,"heatmapProteinLevel.pdf"), width = 10, height = 10)
ggg <- LFQService::plot_heatmap(protintensity("unnest")$data, protintensity("unnest")$config,color = colorRampPalette(c("blue", "white", "red"))(100),
                                breaks = seq(from = -4, to = 4, length.out = 101))
print(ggg)
dev.off()


if (TRUE) {
    library(RColorBrewer)
    data <- protintensity("unnest")$data
    config <- protintensity("unnest")$config
    res <-  toWideConfig(data, config , as.matrix = TRUE)
    
    
    
    newnom <- readxl::read_xlsx("Patient ID Nomenclature for heatmaps.xlsx")
    colnames(newnom) <- make.names(colnames(newnom))
    newnom <- newnom %>% select(1:2)
    colnames(newnom)
    
    annotation_col <- res$annotation
    dim(annotation_col)
    annotation_col <- inner_join(annotation_col, newnom , by = c("Patient_ID" = "Patient.ID.on.proteomics.heatmap"))
    
    annotation_col$Patient_ID <- NULL
    annotation_col <- annotation_col %>% rename(Patient_ID = Patient.ID..final.nomenclature )
    
    Key_for_patient_IDs <- read_excel("data/Key for patient IDs.xlsx")
    names(Key_for_patient_IDs) <- c("Patient_ID", "PatNEW")
    Key_for_patient_IDs <- na.omit(Key_for_patient_IDs)
    Key_for_patient_IDs <- Key_for_patient_IDs %>% distinct()
    
    annotation_colB <- inner_join(annotation_col, Key_for_patient_IDs , by ="Patient_ID")
    dim(annotation_colB)
    annotation_col <- dplyr::rename(annotation_colB, PatOld = Patient_ID, Patient_ID= PatNEW )
    
    annotation_col <- annotation_col %>% mutate(celltype_ = gsub("CMP/MEP", "CMP_MEP", celltype_))
    annotation_col <- rename(annotation_col,
                             `cell type` = celltype_,
                             `cell type` = celltype_,
                             `class therapy progressed` = class_therapy_progressed_ , `patient ID` = Patient_ID)
    resdata <- res$data
    
    #annotation <- lapply(annotation_col , as.character)
    #annotation_col <- data.frame(annotation, stringsAsFactors = FALSE)
    annotation_col[] <- lapply(annotation_col, as.character)
    annotation_col <- annotation_col %>% arrange(`class therapy progressed`,`cell type`,`patient ID`)
  
    ann_colors <- list(
      "cell type" = setNames(brewer.pal(length(unique( annotation_col$`cell type`)), "Set3"), unique( annotation_col$`cell type`)),
      "class therapy progressed" = setNames(brewer.pal(length(unique( annotation_col$`class therapy progressed`)), "Set1"), unique( annotation_col$`class therapy progressed`)),
      "patient ID" = setNames(metafolio::gg_color_hue(length(unique( annotation_col$`patient ID`) ) ), sort(unique( annotation_col$`patient ID`)))
    )
    ann_colors <- readRDS("data/ann_colors222.rds")
    #foo
    
    #ann_colors$`cell type`[c("CMP", "CMP_MEP", "MEP")] <- brewer.pal(9, "YlGn")[4:6]
    annot <- annotation_col
    ann_colors$`cell type` <- ann_colors$`cell type`[-3]
    
    factors <- dplyr::select_at(annot , c('cell type', 'class therapy progressed', 'patient ID'))
    factors <- as.data.frame(factors)
    rownames(factors) <- annot$sampleName
  
    
    resdata <- res$data
    resdata <- quantable::removeNArows(resdata,floor(ncol(resdata)*0.4))
    
    
    #graphics.off()
    # not showing row dendrogram trick
    res <- pheatmap::pheatmap(resdata,
                              scale = "row", silent = TRUE)
    
    res <- pheatmap::pheatmap(resdata[res$tree_row$order,],
                              cluster_rows  = FALSE,
                              scale = "row",
                              annotation_col = factors,
                              show_rownames = F,
                              border_color = NA,
                              annotation_colors = ann_colors,
                              silent = TRUE,
                              color = colorRampPalette(c("blue", "white", "red"))(100),
                              breaks = seq(from = -4, to = 4, length.out = 101))
    
    pdf(file.path(path_qc,"heatmapProteinLevel_222.pdf"), width = 10, height = 10)
    print(res)
    graphics.off()
    
    resdata_CMPMEP <- resdata[,grepl("CMP/MEP",colnames(resdata))]  
    
    res <- pheatmap::pheatmap(resdata_CMPMEP,
                              scale = "row", silent = TRUE)
    
    res <- pheatmap::pheatmap(resdata_CMPMEP[res$tree_row$order,],
                              cluster_rows  = FALSE,
                              scale = "row",
                              annotation_col = factors,
                              show_rownames = F,
                              border_color = NA,
                              annotation_colors = ann_colors,
                              silent = TRUE,
                              color = colorRampPalette(c("blue", "white", "red"))(100),
                              breaks = seq(from = -4, to = 4, length.out = 101))
    pdf(file.path(path_qc,"heatmapProteinLevel_CMP_MEP_222.pdf"), width = 10, height = 10)
    print(res)
    graphics.off()
    
    resdata_HSC <- resdata[,grepl("HSC",colnames(resdata))]  
    
    res <- pheatmap::pheatmap(resdata_HSC,
                              scale = "row", silent = TRUE)
    
    res <- pheatmap::pheatmap(resdata_HSC[res$tree_row$order,],
                              cluster_rows  = FALSE,
                              scale = "row",
                              annotation_col = factors,
                              show_rownames = F,
                              border_color = NA,
                              annotation_colors = ann_colors,
                              silent = TRUE,
                              color = colorRampPalette(c("blue", "white", "red"))(100),
                              breaks = seq(from = -4, to = 4, length.out = 101))
    pdf(file.path(path_qc,"heatmapProteinLevel_HSC_222.pdf"), width = 10, height = 10)
    print(res)
    graphics.off()
  
}

library(ggfortify)

ff <- as.matrix(protintensity("wide")$data[,-1])

rownames(ff) <- unlist( protintensity("wide")$data[,1])
ff <- na.omit(ff)
ff <- t(ff)
ids <- protintensity("wide")$annotation
ids$class_therapy_progressed_ <- gsub("\\.no$|\\.yes$","",ids$class_therapy_progressed_)
#dev.off()


pdf(file.path(path_qc,"PCAproteinLevel.pdf"), width = 6, height = 5)
ids <- ids %>% mutate(celltype_ = gsub("CMP/MEP","CMP_MEP", celltype_))
rownames(ff)
p <- autoplot(prcomp(ff),
         data = ids,
         colour = "celltype_",
         shape = "class_therapy_progressed_", size=5)
p <- p + scale_color_manual(values = ann_colors$`cell type`)
p + theme_light()
dev.off()


xx <- protintensity("plot")
#print(xx$plot[[1]])
#pdf(file.path(path_qc,"AllProteinPlots_WithMissing.pdf"), width = 10, height = 10)
#lapply(xx$plot, print)
#dev.off()


