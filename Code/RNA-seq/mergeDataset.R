# merge o5299
library(ezRun)
library(readr)
library(tibble)
setwd("/scratch/gtan/p2378-Fabienne")

## NovaSeq_20190319_NOV84_o5299_part
dataset1Fn <- "/srv/gstore/projects/p2378/NovaSeq_20190319_NOV84_o5299_part/dataset.tsv"
dataset1 <- read_tsv(dataset1Fn)

## NovaSeq_20190321_NOV85_o5299
dataset2Fn <- "/srv/gstore/projects/p2378/NovaSeq_20190321_NOV85_o5299/dataset.tsv"
dataset2 <- read_tsv(dataset2Fn)

# debug(ezCombineReadDatasets)
datasetCombined <- ezCombineReadDatasets(ezRead.table(dataset1Fn), 
                                         ezRead.table(dataset2Fn),
                                         newDsDir="p2378/NovaSeq_o5299")
write_tsv(as_tibble(datasetCombined, rownames="Name"),
          path=file.path("NovaSeq_o5299", "dataset.tsv"))
