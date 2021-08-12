# annotate and filter vcf file from GTAK call
library(VariantAnnotation)
library(readr)
library(stringr)
library(dplyr)
setwd("/home/gtan/analysis/p2378-Fabienne/RNAVariants_2020-02-16")
vcfFn <- "/srv/gstore/projects/p2378/GatkRnaHaplotyper_37792_2019-07-04--12-31-38/GATK_RnaVariants-haplo.vcf.gz"

## subset vcf for known genes
genes <- c("JAK2", "CALR", "TET2", "ASXL1", "DNMT3A", "MPL", "EZH2", "PPM1D", "NFE2", "SF3B1", "NF1", "TP53", "SRSF2", "U2AF1", "CBL", "KMT2C", "ZRSR2", "GNAS", "IDH2", "SH2B3", "KRAS", "RB1", "PTPN11", "SETBP1", "KIT", "BCOR", "NRAS", "IDH1", "STAG2", "CUX1", "RUNX1", "PHF6", "FLT3", "GATA2", "MBD1", "GNB1")
cosIds <- c("COSM13024","COSM28043","COSM5945302","COSM28040","COSM7344869",
            "COSM6907715","COSM87000","COSM249799","COSM131557","COSM97131",
            "COSM5762885","COSM28026","COSM110752","COSM1594211","COSM1737944",
            "COSM4766107","COSM110785","COSM5020013","COSM3760322","COSM3631397",
            "COSM3762469","COSM1673705","COSM6936295","COSM12600","COSM1651089",
            "COSM7350265","COSM1465373")

anno <- read_tsv("/srv/GT/reference/Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26/Genes/genes_annotation_byGene.txt")
anno <- makeGRangesFromDataFrame(anno, keep.extra.columns = TRUE)

knownGR <- anno[anno$gene_name %in% genes]

tab <- TabixFile(vcfFn)

vcf <- readVcf(tab, "GRCh38", param=knownGR)
head(vcf)
samples(header(vcf))
writeVcf(vcf, "knownGenes.vcf")

## VEP annotated known genes subset
vcf <- readVcf("knownGenes_vep.vcf")
### Keep variants with Cosmic annotation
vcf <- vcf[grepl("COSM", unstrsplit(info(vcf)[,"CSQ"])) & !grepl("LOW", unstrsplit(info(vcf)[,"CSQ"]))]
writeVcf(vcf, filename="knownGenes_vep_COSM.vcf")
GT <- geno(vcf)$GT

rowGenes <- sapply(info(vcf)$CSQ, function(x){paste(stringi::stri_remove_empty(unique(sapply(strsplit(x, "\\|"), "[", 4))), collapse="_")})
COSMIds <- sapply(info(vcf)$CSQ, function(x){paste(unique(unlist(str_extract_all(x, "COSM\\d+"))), collapse=", ")})
loci <- rownames(vcf)
rownames(GT) <- paste(rowGenes, COSMIds)

anno <- read_tsv("FeatureCounts_o5299.tsv")
stopifnot(identical(anno$Name, colnames(GT)))
## fix names
fixNameDT <- read_tsv("/home/gtan/analysis/p2378-Fabienne/Heatmap/fixNames.tsv")
foo <- left_join(anno,
                 fixNameDT, by=c("PatientID [Factor]"="Patient ID on transcriptomics heatmap"))
annotation_col <- data.frame(row.names = foo$Name,
                             "cell type"=sub(".*(NO|HU)_", "", foo$`Condition [Factor]`),
                             "class therapy progressed"=paste(recode(sub("_.*$", "", foo$`Condition [Factor]`), c="control", p="patient"),
                                                              sub("_.*$", "", sub("(p|c)_", "", foo$`Condition [Factor]`)),
                                                              foo$Progressed, sep="."),
                             "patient ID"=foo$`Patient ID, final nomenclature`,
                             check.names = FALSE)
annotation_col$name <- paste(annotation_col$`cell type`,
                             annotation_col$`class therapy progressed`,
                             annotation_col$`patient ID`, sep="~")
colnames(GT) <- annotation_col[colnames(GT), "name"]

library(ComplexHeatmap)
library(RColorBrewer)
textCols <- rep("black", nrow(GT))
textCols[grep("(COSM28043|COSM5945302|COSM28040|COSM1594211|COSM3760322|COSM1673705)", rownames(GT))] <- "blue"
textCols[grep("(COSM1737944|COSM4766107)", rownames(GT))] <- "grey"
ha <- rowAnnotation(foo = anno_text(rownames(GT), gp =gpar(col=textCols)))
p <- Heatmap(GT, col=setNames(RColorBrewer::brewer.pal(length(names(table(GT))), "Set1"), names(table(GT))),
             column_split = annotation_col[match(colnames(GT), annotation_col$name), "patient ID"], name="Genotype",
             row_names_max_width=unit(20, "cm"), right_annotation = ha, show_row_names=FALSE)
p_gg <- grid.grabExpr(draw(p))
cowplot::save_plot(filename="knownGenes.pdf", p_gg, base_height=20, base_width = 40)
### Without Control
patientTokeep <- grep("patient", colnames(GT))
p <- Heatmap(GT[ ,patientTokeep], col=setNames(RColorBrewer::brewer.pal(length(names(table(GT[ ,patientTokeep]))), 
                                                                        "Set1"), names(table(GT[ ,patientTokeep]))),
             column_split = annotation_col[match(colnames(GT), annotation_col$name), "patient ID"][patientTokeep], name="Genotype",
             row_names_max_width=unit(20, "cm"), right_annotation = ha, show_row_names=FALSE)
p_gg <- grid.grabExpr(draw(p))
cowplot::save_plot(filename="knownGenes_Patients.pdf", p_gg, base_height=20, base_width = 30)

## Curated version
loci2Keep <- sapply(relist(unlist(strsplit(COSMIds, ", ")) %in% cosIds,
                           strsplit(COSMIds, ",")), any)
vcf_curated <- vcf[loci2Keep, ]
writeVcf(vcf_curated, filename="knownGenes_vep_COSM_curated.vcf")

GT_curated <- GT[loci2Keep, ]
textCols <- rep("black", nrow(GT_curated))
textCols[grep("(COSM28043|COSM5945302|COSM28040|COSM1594211|COSM3760322|COSM1673705)", rownames(GT_curated))] <- "blue"
textCols[grep("(COSM1737944|COSM4766107)", rownames(GT_curated))] <- "grey"
ha <- rowAnnotation(foo = anno_text(rownames(GT_curated), gp =gpar(col=textCols)))
p <- Heatmap(GT_curated, col=setNames(RColorBrewer::brewer.pal(length(names(table(GT))), "Set1"), names(table(GT))),
             column_split = annotation_col[match(colnames(GT_curated), annotation_col$name), "patient ID"], name="Genotype",
             row_names_max_width=unit(20, "cm"), right_annotation = ha, show_row_names=FALSE)
p_gg <- grid.grabExpr(draw(p))
cowplot::save_plot(filename="knownGenes_curated.pdf", p_gg,
                   base_height=10, base_width = 40)
### Without Control
patientTokeep <- grep("patient", colnames(GT_curated))
p <- Heatmap(GT_curated[ ,patientTokeep], col=setNames(RColorBrewer::brewer.pal(length(names(table(GT_curated[ ,patientTokeep]))), 
                                                                        "Set1"), names(table(GT_curated[ ,patientTokeep]))),
             column_split = annotation_col[match(colnames(GT_curated), annotation_col$name), "patient ID"][patientTokeep], name="Genotype",
             row_names_max_width=unit(20, "cm"), right_annotation = ha, show_row_names=FALSE)
p_gg <- grid.grabExpr(draw(p))
cowplot::save_plot(filename="knownGenes_curated_Patients.pdf", p_gg, base_height=10, base_width = 30)

## Alle fraction
toPlot <- matrix(sapply(geno(vcf_curated)$AD, "[", 2), nrow=nrow(geno(vcf_curated)$AD),
                 ncol=ncol(geno(vcf_curated)$AD), dimnames=dimnames(geno(vcf_curated)$AD)) / 
  geno(vcf_curated)$DP
rownames(toPlot) <- rownames(GT_curated)
colnames(toPlot) <- colnames(GT_curated)

textCols <- rep("black", nrow(GT_curated))
textCols[grep("(COSM28043|COSM5945302|COSM28040|COSM1594211|COSM3760322|COSM1673705)", rownames(GT_curated))] <- "blue"
textCols[grep("(COSM1737944|COSM4766107)", rownames(GT_curated))] <- "grey"
ha <- rowAnnotation(foo = anno_text(rownames(GT_curated), gp =gpar(col=textCols)))

p <- Heatmap(toPlot, column_split = annotation_col[match(colnames(GT_curated), annotation_col$name), "patient ID"], 
             name="Allel Fraction",
             row_names_max_width=unit(20, "cm"), cluster_rows=FALSE,
             cluster_columns=FALSE, right_annotation = ha, show_row_names=FALSE)
p_gg <- grid.grabExpr(draw(p))
cowplot::save_plot(filename="knownGenes_AllelFraction_curated.pdf", p_gg,
                   base_height=10, base_width = 40)
### Without Control
patientTokeep <- grep("patient", colnames(toPlot))
p <- Heatmap(toPlot[, patientTokeep], column_split = annotation_col[match(colnames(GT_curated), annotation_col$name), "patient ID"][patientTokeep], 
             name="Allel Fraction",
             row_names_max_width=unit(20, "cm"), cluster_rows=FALSE,
             cluster_columns=FALSE, right_annotation = ha, show_row_names=FALSE)
p_gg <- grid.grabExpr(draw(p))
cowplot::save_plot(filename="knownGenes_AllelFraction_curated_Patients.pdf", p_gg,
                   base_height=10, base_width = 30)
