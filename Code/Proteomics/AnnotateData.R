rm(list = ls())
library(tidyverse)


tmp5 <- read_tsv("data/20180923_230031_SA_1807FMA_HSPCPV5frwoisolw_SNP_Report.xls")
grepl("EG.PrecursorId", colnames(tmp5))



annotation <- readxl::read_xlsx("Annotation/Proteomics Samples for Data Analysis_Annotation_final_without_1samples.xlsx")
colnames(annotation) <- tolower(make.names(colnames(annotation)))

with(annotation, table(patient.age))
with(annotation, table(diagnosis))
with(annotation, table(celltype))
with(annotation, table(class , diagnosis))
with(annotation, table(class , therapy))

annotation$class_therapy <- with(annotation, interaction(class , therapy))
annotation$class_therapy_progressed <- with(annotation, interaction(class , therapy, progressed))

mean(annotation$filename %in% tmp5$R.FileName)

dataSet <- inner_join(annotation, tmp5, by = c("filename" = "R.FileName"))
nrow(tmp5) == nrow(dataSet)
annotation$filename %in% unique(tmp5$R.FileName)
setdiff(unique(tmp5$R.FileName),annotation$filename)

#resPepProtAnnot <- read_rds(path = "data/20180923_230031_SA_1807FMA_HSPCPV5frwoisolw_SNP_Report.rds")

dataSet$isotope <- "light"
dataSet <- dataSet %>% dplyr::filter(celltype == "HSC" | celltype == "CMP/MEP")


write_rds(dataSet, path = "data/20180923_230031_SA_1807FMA_HSPCPV5frwoisolw_SNP_Report.rds")


