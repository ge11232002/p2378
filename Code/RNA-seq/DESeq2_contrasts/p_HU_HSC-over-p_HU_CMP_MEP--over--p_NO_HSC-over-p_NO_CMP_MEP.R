## The comparison of difference: p_HU_HSC-over-p_HU_CMP_MEP--over--p_NO_HSC-over-p_NO_CMP_MEP
setwd("/home/gtan/analysis/p2378-Fabienne/DESeq2_contrasts")
library(ezRun)
library(stringr)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = 'fgcz-c-042,fgcz-c-045,fgcz-c-046,fgcz-c-048,fgcz-c-056,fgcz-c-063,fgcz-c-065,fgcz-h-004,fgcz-h-006,fgcz-h-007,fgcz-h-008,fgcz-h-009,fgcz-h-010,fgcz-h-011,fgcz-h-012,fgcz-h-013,fgcz-h-014,fgcz-h-015,fgcz-h-016,fgcz-h-017,fgcz-h-018,fgcz-h-019'
param[['process_mode']] = 'DATASET'
param[['samples']] = 'Haemo_CMP,Haemo_CMP_MEP,Haemo_GMP,Haemo_HSC,Haemo_MEP,PID106_CMP_MEP,PID106_GMP,PID106_HSC,PID115_CMP,PID115_CMP_MEP,PID115_GMP,PID115_HSC,PID115_MEP,PID115star_CMP,PID115star_CMP_MEP,PID115star_GMP,PID115star_HSC,PID115star_MEP,PID153_CMP,PID153_CMP_MEP,PID153_GMP,PID153_HSC,PID153_MEP,PID199_GMP,PID206_CMP_MEP,PID206_GMP,PID206_HSC,PID240_GMP,PID282_CMP,PID282_CMP_MEP,PID282_GMP,PID282_HSC,PID282_MEP,PID291_CMP,PID291_CMP_MEP,PID291_GMP,PID291_HSC,PID291_MEP,PID304_CMP,PID304_CMP_MEP,PID304_GMP,PID304_HSC,PID304_MEP,PID318_CMP_MEP,PID318_GMP,PID318_HSC,PID473_CMP,PID473_CMP_MEP,PID473_GMP,PID473_HSC,PID473_MEP,PID480_CMP,PID480_CMP_MEP,PID480_GMP,PID480_HSC,PID480_MEP,PID64_GMP,PID730_CMP,PID730_CMP_MEP,PID730_GMP,PID730_HSC,PID730_MEP,PID757_CMP_MEP,PID757_GMP,PID771_CMP,PID771_CMP_MEP,PID771_GMP,PID771_HSC,PID771_MEP,PID781_CMP,PID781_CMP_MEP,PID781_GMP,PID781_HSC,PID781_MEP,PID81_CMP_MEP,PID81_GMP,PID81_HSC,PID872_CMP,PID872_CMP_MEP,PID872_GMP,PID872_HSC,PID872_MEP,PID92_CMP,PID92_CMP_MEP,PID92_GMP,PID92_HSC,PID92_MEP,PID94_CMP_MEP,PID94_GMP,PID94_HSC'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['grouping']] = 'Condition'
param[['sampleGroup']] = 'p_HU_HSC'
param[['refGroup']] = 'p_NO_HSC'
param[['onlyCompGroupsHeatmap']] = 'false'
param[['runGO']] = 'true'
param[['grouping2']] = 'PatientID'
param[['backgroundExpression']] = '10'
param[['transcriptTypes']] = 'protein_coding'
param[['specialOptions']] = 'doPrecomputeEnrichr=false'
param[['expressionName']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['Rversion']] = 'Dev/R/3.5.1'
param[['comparison']] = 'p_HU_HSC-over-p_HU_CMP_MEP--over--p_NO_HSC-over-p_NO_CMP_MEP' ## ToAdapt
param[['name']] = 'p_HU_HSC-over-p_HU_CMP_MEP--over--p_NO_HSC-over-p_NO_CMP_MEP'## ToAdapt
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2378/DESeq2_35913_p_HU_MEP--over--p_HU_CMP_2019-04-23--14-21-06'
output = list()
output[['Name']] = 'p_HU_HSC-over-p_HU_CMP_MEP--over--p_NO_HSC-over-p_NO_CMP_MEP'## ToAdapt
output[['Species']] = 'Homo sapiens (human)'
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_91-2018-02-26'
output[['Static Report [Link]']] = 'p2378/DESeq2_35913_p_HU_MEP--over--p_HU_CMP_2019-04-23--14-21-06/p_HU_MEP--over--p_HU_CMP/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreDEG_app/?data=p2378/DESeq2_35913_p_HU_MEP--over--p_HU_CMP_2019-04-23--14-21-06/p_HU_MEP--over--p_HU_CMP/result-p_HU_MEP--over--p_HU_CMP-xfkaubexfqmi-EzResult.RData'
output[['Report [File]']] = 'p2378/DESeq2_35913_p_HU_MEP--over--p_HU_CMP_2019-04-23--14-21-06/p_HU_MEP--over--p_HU_CMP'
input = '/srv/gstore/projects/p2378/DESeq2_35913_p_HU_MEP--over--p_HU_CMP_2019-04-23--14-21-06/input_dataset.tsv'
param = ezParam(param)
input = EzDataset$new(file=input, dataRoot=param$dataRoot)
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
# EzAppDeseq2$new()$run(input=input, output=output, param=param)

cwd <- getwd()
setwdNew(param[['name']])
rawData = loadCountDataset(input, param)
param$grouping = input$getColumn(param$grouping)
param$grouping2 = input$getColumn(param$grouping2)
param$sampleGroup <- c("p_HU_HSC", "p_HU_CMP_MEP") ## ToAdapt
param$refGroup <- c("p_NO_HSC", "p_NO_CMP_MEP") ## ToAdapt

## twoGroupCountComparison
x = assays(rawData)$counts
presentFlag = assays(rawData)$presentFlag
job = ezJobStart("twoGroupCountComparison")
metadata(rawData)$analysis <- "NGS two group analysis"
param$testMethod = "deseq2"

metadata(rawData)$method <- param$testMethod

isSample = param$grouping %in% param$sampleGroup
isRef = param$grouping %in% param$refGroup
## compute which probes are present
isPresent = ezPresentFlags(x, presentFlag=presentFlag, param=param,
                           isLog=FALSE)
useProbe = logical(nrow(x))
useProbe[rowMeans(isPresent[, isRef, drop=FALSE]) >= 0.5] = TRUE
useProbe[rowMeans(isPresent[, isSample, drop=FALSE]) >= 0.5] = TRUE
rowData(rawData)$isPresentProbe <- useProbe
assays(rawData)$isPresent <- isPresent
isPresent=useProbe
## run deseq2
x <- round(x)
grouping <- param$grouping
subject <- colData(rawData)$`PatientID [Factor]`
sampleGroup <- param$sampleGroup
refGroup <- param$refGroup
grouping2=param$grouping2

require(DESeq2)
## get size factors -- grouping2 not needed
colData = data.frame(grouping=as.factor(grouping), 
                     row.names=colnames(x))
dds = DESeqDataSetFromMatrix(countData=x, colData=colData,
                             design= ~ grouping)
dds = estimateSizeFactors(dds, controlGenes=isPresent)
sf = 1/dds@colData$sizeFactor

## remove the samples that do not participate in the comparison
isSample = grouping %in% sampleGroup
isRef = grouping %in% refGroup
grouping = grouping[isSample|isRef]
x2 = x[ ,isSample|isRef]
if (ezIsSpecified(grouping2)){
  grouping2 = grouping2[isSample|isRef]
}

colData = data.frame(grp=str_extract(grouping, "(p|c)_(HU|NO)"), ## ToAdapt
                     cnd=sub("(p|c)_(HU|NO)_", "", grouping), ## ToAdapt
                     row.names=colnames(x2))
PatientFactor <- lapply(split(grouping2, colData$grp),
                        as.factor)
PatientFactor <- lapply(PatientFactor, function(x){levels(x) <- 1:length(levels(x));x})
PatientFactor <- unlist(PatientFactor)
names(PatientFactor) <- sub("(p|c)_(HU|NO).", "", names(PatientFactor))
colData$ind.n <- PatientFactor[rownames(colData)]
colData <- colData[!colData$ind.n %in% c("8"), ] ## ToAdapt otherwise it complains not full rank
colData$ind.n <- factor(colData$ind.n)

dds = DESeqDataSetFromMatrix(countData=x2[ ,rownames(colData)],
                             colData=colData,
                             design= ~ cnd + grp)

design(dds) <- ~ grp + grp:ind.n + grp:cnd
dds = estimateSizeFactors(dds, controlGenes=isPresent)
dds = DESeq(dds)
resultsNames(dds)
res = results(dds, contrast=list("grpp_HU.cndHSC", "grpp_NO.cndHSC"))
res = as.list(res)
res$sf = sf

## twoGroupCountComparison
rowData(rawData)$log2Ratio <- res$log2FoldChange
metadata(rawData)$fitGlm = res$fitGlm
colData(rawData)$sf <- res$sf
pValue = res$pval
pValue[is.na(pValue)] = 1
metadata(rawData)$nativeResult <- res
useProbe[is.na(useProbe)] = FALSE
fdr = rep(NA, length(pValue))
fdr[useProbe] = p.adjust(pValue[useProbe], method="fdr")
rowData(rawData)$pValue <- pValue
rowData(rawData)$fdr <- fdr
## usedInTest was used before to pre-filter. Not used for now.
rowData(rawData)$usedInTest = useProbe
assays(rawData)$xNorm = ezScaleColumns(x, colData(rawData)$sf)

metadata(rawData)$summary = c("Name"=param$name,
                              "Reference Build"=param$refBuild,
                              "Feature Level"=metadata(rawData)$featureLevel,
                              "Normalization"=param$normMethod)
param$refGroup <- param$refGroup[1]
param$sampleGroup <- param$sampleGroup[1]
metadata(rawData)$param <- param
deResult = EzResult(se=rawData)

styleFiles <- file.path(system.file("templates", package="ezRun"),
                        c("fgcz.css", "twoGroups.Rmd",
                          "fgcz_header.html", "banner.png"))
file.copy(from=styleFiles, to=".", overwrite=TRUE)
rmarkdown::render(input="twoGroups.Rmd", envir = new.env(),
                  output_dir=".", output_file="00index.html", quiet=TRUE)

