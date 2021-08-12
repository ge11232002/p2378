# heatmap from CountQC
library(SummarizedExperiment)
library(ezRun)
library(readr)
library(dplyr)
library(pheatmap)
library(cowplot)
setwd("/home/gtan/analysis/p2378-Fabienne/Heatmap")

## fix names
fixNameDT <- read_tsv("fixNames.tsv")
fixNameDT2 <- read_tsv("fixName2.tsv")
fixNameDT <- left_join(fixNameDT, fixNameDT2, by=c("Patient ID, final nomenclature"="Old ID"))

## CountQC
# https://fgcz-sushi.uzh.ch/data_set/p2378/35914
load("/srv/gstore/projects/p2378/CountQC_35913_2019-04-04--12-12-55/Count_QC/counts-olbcexapyxpk-EzResult.RData")
rawData <- se

output <- metadata(rawData)$output
param <- metadata(rawData)$param
seqAnno <- data.frame(rowData(rawData), row.names=rownames(rawData),
                      check.names = FALSE, stringsAsFactors=FALSE)
dataset <- data.frame(colData(rawData), 
                      row.names=colnames(rawData), check.names = FALSE,
                      stringsAsFactors=FALSE)


types = data.frame(row.names=rownames(seqAnno))
for (nm in setdiff(na.omit(unique(seqAnno$type)), "")){
  types[[nm]] = seqAnno$type == nm
}

design = ezDesignFromDataset(dataset, param)
samples = rownames(design)
nSamples = length(samples)
conds = ezConditionsFromDesign(design, maxFactors = 2)
nConds = length(unique(conds))
sampleColors = getSampleColors(conds)

signal = shiftZeros(getSignal(rawData), param$minSignal)
presentFlag = assays(rawData)$presentFlag
signalRange = range(signal, na.rm=TRUE)
log2Signal = log2(signal)
isPresent = ezPresentFlags(assays(rawData)$counts, presentFlag=presentFlag, 
                           param=param, isLog=metadata(rawData)$isLog)
signalCond = 2^averageColumns(log2Signal, by=conds)
isPresentCond = averageColumns(isPresent, by=conds) >= 0.5
isPresentStudy <- rowMeans(isPresentCond) >= 0.5

assays(rawData)$signal = signal

isValid = isPresentStudy
if (!is.null(seqAnno$IsControl)){
  isValid = isValid & !seqAnno$IsControl
}

if(sum(isValid) < 10){
  cat("Not enough valid features for further plots", "\n")
  knit_exit()
}

x = log2(2^log2Signal[isValid, ] + param$backgroundExpression)
xNormed = sweep(x, 1 , rowMeans(x))
xSd <- setNames(rowSds(x, na.rm=TRUE), rownames(x))
ord = order(xSd, decreasing=TRUE)
topGenes = ord[1:min(param$topGeneSize, length(ord))]

sdThresh = param$highVarThreshold
use = xSd > sdThresh & !rowAnyMissings(x)

if (sum(use, na.rm=TRUE) > param$maxGenesForClustering){
  use[use] = rank(-1 * xSd[use], ties.method="max") <= param$maxGenesForClustering
  sdThresh = signif(min(xSd[use]), digits=3)
}

## pheatmap setup
callback = function(hc, mat){
  dend = reorder(as.dendrogram(hc), wts = rowMeans(mat, na.rm=TRUE))
  as.hclust(dend)
}
lim=c(-4, 4)
colors=getBlueRedScale()

## Including BC samples
toPlot <- xNormed[use, ]
foo <- left_join(as_tibble(design, rownames = "id"),
                 fixNameDT, by=c("PatientID"="Patient ID on transcriptomics heatmap"))
annotation_col <- data.frame(row.names = foo$id,
                             "cell type"=sub(".*(NO|HU)_", "", foo$Condition),                             
                             "class therapy progressed"=paste(recode(sub("_.*$", "", foo$Condition), c="control", p="patient"),
                                                              sub("_.*$", "", sub("(p|c)_", "", foo$Condition)),
                                                              foo$Progressed, sep="."),
                             "patient ID"=foo$`New ID`,
                             check.names = FALSE)
annotation_col$name <- paste(annotation_col$`cell type`,
                             annotation_col$`class therapy progressed`,
                             annotation_col$`patient ID`, sep="~")

ann_colors <- readRDS("ann_colors.rds")
ann_colors$`patient ID` <- ann_colors$`patient ID`[names(ann_colors$`patient ID`) %in% annotation_col$`patient ID`]

colnames(toPlot) <- annotation_col[colnames(toPlot), "name"]
rownames(annotation_col) <- annotation_col$name
annotation_col$name <- NULL

pheatmap(toPlot, color=colors, clustering_method="ward.D2",
         breaks=seq(from=lim[1], to=lim[2], length.out=257),
         scale="none", show_rownames = FALSE,
         clustering_callback = callback,
         treeheight_row=0, annotation_col=annotation_col,
         annotation_colors = ann_colors,
         filename="cluster-heatmap_withoutBC.pdf",
         height=10, width=15
         )

## MDS plot
design2 <- data.frame(row.names = rownames(design),
                      "cell type"=sub(".*(NO|HU)_", "", design$Condition),
                      "class therapy"=paste(recode(sub("_.*$", "", design$Condition), c="control", p="patient"),
                                            sub("_.*$", "", sub("(p|c)_", "", design$Condition)),  sep="."),
                      check.names = FALSE
)
ezMdsGG2 <- function(signal, design, ndim=2, main="MDS plot"){
  require(edgeR)
  require(plotly)
  require(ggrepel)
  if(ndim != 2){
    stop("ggplot2 only produces 2D plot")
  }
  
  y = DGEList(counts=signal, group=colnames(signal))
  mds = plotMDS(y, plot=FALSE, ndim=ndim)
  toPlot <- data.frame(samples=colnames(signal),
                       design,
                       stringsAsFactors = FALSE,
                       check.names = FALSE)
  mdsOut <- mds$cmdscale.out
  colnames(mdsOut) <- c("Leading logFC dim1", "Leading logFC dim2")
  toPlot <- cbind(toPlot, mdsOut)
  if(ncol(design) > 1L){
    p <- ggplot(toPlot, aes(`Leading logFC dim1`, `Leading logFC dim2`)) +
      geom_point(aes_string(colour=sym(colnames(design)[1]),
                            shape =sym(colnames(design)[2])),
                 size = 3) +
      ggtitle(main)
  }else{
    p <- ggplot(toPlot, aes(`Leading logFC dim1`, `Leading logFC dim2`)) +
      geom_point(aes_string(colour=sym(colnames(design)[1])),
                 size = 3) + 
      ggtitle(main)
  }
  p
}
p <- ezMdsGG2(signal=x, design=design2, ndim=2)
p <- p + scale_colour_manual(values=ann_colors$`cell type`) +
  theme_half_open() + background_grid()
save_plot(filename="MDS_PresentGenes_withoutBC.pdf", p, base_height = 5)

save.image("withoutBC.RData")
