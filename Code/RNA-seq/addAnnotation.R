library(readr)
library(readxl)
library(dplyr)

setwd("/Users/gtan/OneDrive/Projects-FGCZ/FGCZ/p2378-Fabienne/Notes")

## Add anno to FeatureCounts
featCounts <- read_tsv("FeatureCounts_35699.tsv")
featCounts2 <- read_tsv("FeatureCounts_35724.tsv")
featCounts <- filter(featCounts, Name !="PID291_GMP") %>%
  bind_rows(featCounts2) %>%
  arrange(Name)
anno <- read_excel("RNAseq Annotation File.xlsx")

featCounts <- left_join(featCounts, anno, by=c("Name"="Sample Name"))
featCounts <- select(featCounts, -Condition) %>%
  mutate(Celltype=sub("/", "_", Celltype)) %>%
  mutate(`Condition [Factor]`=paste(Class, Treatment, Celltype, sep="_")) %>%
  mutate(`Patient ID`=make.names(`Patient ID`)) %>%
  rename("PatientID [Factor]"=`Patient ID`)
featCounts <- select(featCounts, Name, `Condition [Factor]`,
                     `PatientID [Factor]`, Class, Treatment, Celltype,
                     Diagnosis, everything())
write_tsv(featCounts, path="FeatureCounts_o5299.tsv")

## Add anno to STAR
dataset1 <- read_tsv("STAR_35662.tsv")
dataset2 <- read_tsv("STAR_35717.tsv")
dataset <- filter(dataset1, Name !="PID291_GMP") %>%
  bind_rows(dataset2) %>%
  arrange(Name)
anno <- read_excel("RNAseq Annotation File.xlsx")

dataset <- left_join(dataset, anno, by=c("Name"="Sample Name"))
dataset <- select(dataset, -Condition) %>%
  mutate(Celltype=sub("/", "_", Celltype)) %>%
  mutate(`Condition [Factor]`=paste(Class, Treatment, Celltype, sep="_")) %>%
  mutate(`Patient ID`=make.names(`Patient ID`)) %>%
  rename("PatientID [Factor]"=`Patient ID`)
dataset <- select(dataset, Name, `Condition [Factor]`,
                     `PatientID [Factor]`, Class, Treatment, Celltype,
                     Diagnosis, everything())
write_tsv(dataset, path="STAR_o5299.tsv")

## Add anno to Kallisto
featCounts <- read_tsv("Kallisto_35662.tsv")
featCounts2 <- read_tsv("Kallisto_35717.tsv")
featCounts <- filter(featCounts, Name !="PID291_GMP") %>%
  bind_rows(featCounts2) %>%
  arrange(Name)
anno <- read_excel("RNAseq Annotation File.xlsx")

featCounts <- left_join(featCounts, anno, by=c("Name"="Sample Name"))
featCounts <- select(featCounts, -Condition) %>%
  mutate(Celltype=sub("/", "_", Celltype)) %>%
  mutate(`Condition [Factor]`=paste(Class, Treatment, Celltype, sep="_")) %>%
  mutate(`Patient ID`=make.names(`Patient ID`)) %>%
  rename("PatientID [Factor]"=`Patient ID`)
featCounts <- select(featCounts, Name, `Condition [Factor]`,
                     `PatientID [Factor]`, Class, Treatment, Celltype,
                     Diagnosis, everything())
write_tsv(featCounts, path="Kallisto_o5299.tsv")
