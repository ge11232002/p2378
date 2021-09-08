library(readr)
library(fgcz.gsea.ora)
library(tidyverse)

genomics <- readr::read_tsv("Count_QC-rpkm.txt")
proteomics <- readxl::read_xlsx("second_Draft_V4_WithContaminants/path_qc/transformed_ProteinIntensities_WIDE.xlsx")


colnames(genomics)
colnames(proteomics)
proteomics$row_Id <- 1:nrow(proteomics)
dim(proteomics)
dim(proteomics)
prot2 <- proteomics %>% tidyr::separate_rows(protein_Id, sep=";")

dim(prot2)
sum(duplicated(prot2$row_Id))

#quantable::write.vector(proteomics$protein_Id, file = "tmp.txt")


prot_mapped_table <- fgcz.gsea.ora::map_ids_uniprot(prot2,ID_col = "protein_Id", to = "ENSEMBL_ID")
plot(table(prot_mapped_table$protein_Id))
xx <- sort(as.array(table(prot_mapped_table$protein_Id)), decreasing = TRUE)
(xx[xx > 1])


dim(prot2)
dim(prot_mapped_table)

prot_mapped_table$protein_Id[is.na(prot_mapped_table$ENSEMBL_ID)]



prot_mapped_table$ENSEMBL_ID[!mapped_table$ENSEMBL_ID %in% genomics$`Feature ID`]



sum(mapped_table$ENSEMBL_ID %in% genomics$`Feature ID`)
left_join(mapped_table,genomics, by = c("ENSEMBL_ID" = "Feature ID"))

