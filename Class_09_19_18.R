#class

library(affy)


setwd("/Users/Amanda/Desktop/estrogen")


targetsFile <- "estrogen.txt"


pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)

ER <- pData(pd)$estrogen

Time <- factor(pData(pd)$time.h)


design <- model.matrix(~ER+Time)

# ~ is a formula. Additive effect / no interaction. Goes to one time point to determine significance


design2 <- model.matrix(~ER*Time)

design2

#multiply to see interaction. Goes to both time points to determine significance.

design
design

raw <-ReadAffy(celfile.path = "C:/Users/Bio-user/Documents/GitHub/Class-Exercise-1/estrogen", filenames=rownames(pData(pd)),phenoData = pd)


raw

eset <- rma(raw)

library(limma)

fit1 <- lmFit(eset, design)
# fit of expression values on a linear line


fit1 <- eBayes(fit1)
#e bayes function to calculate p value and fold change


topTable(fit1, coef=2)
#get top probes ids with significant value

fit2 <- lmFit(eset, design2)

fit2 <- eBayes(fit2)

fit2 <- eBayes(fit2)

topTable(fit2, coef=2)

#use both additive and interaction and find common genes

#Annotation of gene names

source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

library(genefilter)


library(GEOquery)


library(limma)


url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz"


filenm <- "GSE33126_series_matrix.txt.gz"


if(!file.exists(filenm)) download.file(url, destfile=filenm)


gse <- getGEO(filename=filenm)

#varfilter, fData, anno

gse.expfilter <- varFilter(gse)

anno <- fData(gse.expfilter)

anno <- anno[,c("Symbol","Entrez_Gene_ID","Chromosome","Cytoband")]

fit2$genes <- anno

topTable(fit2)

#make volcanoplot

volcanoplot(fit2)

