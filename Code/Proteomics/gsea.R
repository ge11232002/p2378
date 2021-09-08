rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(org.Hs.eg.db)
library(conflicted)

fpath <-"interaction_draft/Contrasts_SignificanceValues_f_Cells_Treatment_Interaction_PIVOT.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))
#dd %>% separate(protein_Id, c("sp","uniid","trash"),sep="\\|", remove = FALSE)-> dd
#dd %>% dplyr::filter(!is.na(uniid)) -> dd

colnames(dd)

organism <- "hsapiens"
ID_col <- "protein_Id"

columns <- grep("estimate",colnames(dd), value=TRUE)


for( i in columns){
  print(i)  
  if(!dir.exists(i)){
    dir.create(i)
  }
  
  fc_col <- i
  target <- "geneontology_Biological_Process_noRedundant"
  map_col <- "GO"
  
  
  ranktable <- dplyr::select(dd, ID = !!sym(ID_col), Score = !!sym(fc_col))
  
  GSEA_res <-
    WebGestaltR(
      enrichMethod = "GSEA",
      # does permutation test with 1000 permutations per default, might take a while
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable,
      interestGeneType = "uniprotswissprot",
      # Or "uniprotswissprot" in some cases
      outputDirectory = fc_col,
      isOutput = TRUE,
      perNum = 100,
      projectName = "GSEA_proj"
    )
}

