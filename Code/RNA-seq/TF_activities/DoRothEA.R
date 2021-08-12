library(dorothea)
library(bcellViper)
library(viper)
library(tidyverse)
library(SummarizedExperiment)

setwd("/home/gtan/analysis/p2378-Fabienne/TF_activities")

data(dorothea_hs, package = "dorothea")
write_tsv(dorothea_hs  %>% filter(target == "PF4"),
          "PF4_TFs.tsv")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))
         
## Load the expression data
# https://fgcz-gstore.uzh.ch/projects/p2378/DESeq2_35913_p_NO_HSC--over--c_NO_HSC_2019-04-04--13-11-10/p_NO_HSC--over--c_NO_HSC/00index.html
load("/srv/gstore//projects/p2378/DESeq2_35913_p_NO_HSC--over--c_NO_HSC_2019-04-04--13-11-10/p_NO_HSC--over--c_NO_HSC/result-p_NO_HSC--over--c_NO_HSC-rwzchdpxgvpe-EzResult.RData")
se <- se[ ,se$`Condition [Factor]` %in% c("p_NO_HSC", "c_NO_HSC", "p_HU_HSC")]
se <- se[ ,order(se$`Condition [Factor]`, decreasing = TRUE)]
exprs <- assay(se, "xNorm")
rownames(exprs) <- rowData(se)$gene_name

tf_activities <- run_viper(exprs, dorothea_hs, 
                           options =  list(method = "scale", minsize = 4,
                                           eset.filter = FALSE, cores = 1,
                                           verbose = FALSE))
tf4PF4 <- dorothea_hs  %>% filter(target == "PF4") %>% pull(tf) %>%
  c("GATA1", "RUNX1", "FLI1", "ELF1", "GABPA", "TAL1")
write_tsv(as_tibble(tf_activities[tf4PF4, ], rownames="gene_name"), 
          "tf_activities.tsv")
