
setwd("/home/gtan/analysis/p2378-Fabienne/CountQC")
library(readr)
library(dplyr)

# https://fgcz-sushi.uzh.ch/data_set/p2378/35761
fpkm <- read_tsv("/srv/gstore/projects/p2378/CountQC_35760_2019-03-31--16-13-38/Count_QC/Count_QC-rpkm.txt")

input_dataset <- read_tsv("/srv/gstore/projects/p2378/CountQC_35760_2019-03-31--16-13-38/input_dataset.tsv")

fixNameDT <- read_tsv("/home/gtan/analysis/p2378-Fabienne/Heatmap/fixNames.tsv")

dataset <- left_join(input_dataset, fixNameDT, by=c("PatientID [Factor]"="Patient ID on transcriptomics heatmap"))

annotation_col <- tibble(Name = dataset$Name,
                         "cell type"=sub(".*(NO|HU)_", "", dataset$`Condition [Factor]`),                             
                         "class therapy progressed"=paste(recode(sub("_.*$", "", dataset$`Condition [Factor]`), c="control", p="patient"),
                                                              sub("_.*$", "", sub("(p|c)_", "", dataset$`Condition [Factor]`)),
                                                              dataset$Progressed, sep="."),
                         "patient ID"=dataset$`Patient ID, final nomenclature`)
annotation_col$samename <- paste(annotation_col$`cell type`,
                             annotation_col$`class therapy progressed`,
                             annotation_col$`patient ID`, sep="~")

stopifnot(setdiff(colnames(fpkm), "Feature ID") %in% annotation_col$Name)

mapping <- setNames(annotation_col$Name, annotation_col$samename)

fpkm <- dplyr::rename(fpkm, !!mapping)
write_tsv(fpkm, "Count_QC-rpkm.txt")
