
library(VariantAnnotation)
fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcf <- readVcf(fl, "hg19")
vcf

header(vcf)
samples(header(vcf))
geno(header(vcf))

head(rowRanges(vcf), 3)
ref(vcf)[1:5]
qual(vcf)[1:5]
alt(vcf)[1:5]

geno(vcf)
sapply(geno(vcf), class)

info(vcf)[1:4, 1:5]

library(SNPlocs.Hsapiens.dbSNP151.GRCh38)
rd <- rowRanges(vcf)
ch22snps <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP151.GRCh38, "22")

dbsnpchr22 <- names(rd) %in% ch22snps$RefSNP_id
table(dbsnpchr22)

info(header(vcf))[c("VT", "LDAF", "RSQ"),]

metrics <- data.frame(QUAL=qual(vcf), inDbSNP=dbsnpchr22,
                      VT=info(vcf)$VT, LDAF=info(vcf)$LDAF, 
                      RSQ=info(vcf)$RSQ)

library(ggplot2)
ggplot(metrics, aes(x=RSQ, fill=inDbSNP)) +
       geom_density(alpha=0.5) +
       scale_x_continuous(name="MaCH / Thunder Imputation Quality") +
       scale_y_continuous(name="Density") +
       theme(legend.position="top")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(vcf) <- "chr22"
rd <- rowRanges(vcf)
loc <- locateVariants(rd, txdb, CodingVariants())
head(loc, 3)

