# code pulled from vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start

# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("apeglm")
install.packages("corrplot")

BiocManager::install("gplots")

# ---------------- IMPORTS ---------------- #

# import packages
library("pasilla")
library("DESeq2")
library("apeglm")
library("corrplot")
library("dplyr")
library("RColorBrewer")
library("gplots")

# ---------------- DATA ---------------- #

setwd("~/Documents/Github/DESeq2_Tutorial/")

# get data
cts <- as.matrix(read.csv(file = "data/gc_ancg_b.csv", row.names="refseq"))
coldata <- read.csv("data/ancg_b.csv", row.names="Run")
coldata <- coldata[,c("source_name","age_at_death", "Brain_pH", "post.mortem_interval")]
coldata$source_name <- factor(coldata$source_name)
coldata$age_at_death <- factor(coldata$age_at_death)
coldata$Brain_pH <- factor(coldata$Brain_pH)
coldata$post.mortem_interval <- factor(coldata$post.mortem_interval)
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# ---------------- DESeqDataSet ---------------- #

# create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~source_name)

# set up metadata
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)

# ---------------- FILTERING ---------------- #

# pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# factors in R!
dds$source_name <- factor(dds$source_name, levels = c("AnCg_Control", "AnCg_Bipolar Disorder"))

# ---------------- DEA ---------------- #

# carries out: estimation of size factors, estimation of dispersion: neg. binomial GLM
dds <- DESeq(dds)
#dds
res <- results(dds)

resname <- resultsNames(dds)
resname[2]

# Log fold change
resLFC <- lfcShrink(dds, coef=resname[2], type="apeglm")
resLFC

resOrdered <- resLFC[order(resLFC$pvalue),]
sum(resOrdered$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

res_ancg_b <- res
res_ancg_m <- res
res_ancg_s <- res
res_nacc_b <- res
res_nacc_m <- res
res_nacc_s <- res
res_dlpfc_b <- res
res_dlpfc_m <- res
res_dlpfc_s <- res

write.csv(res_ancg_b, "corrplot/ancg_b.csv")
write.csv(res_ancg_m, "corrplot/ancg_m.csv")
write.csv(res_ancg_s, "corrplot/ancg_s.csv")
write.csv(res_nacc_b, "corrplot/nacc_b.csv")
write.csv(res_nacc_m, "corrplot/nacc_m.csv")
write.csv(res_nacc_s, "corrplot/nacc_s.csv")
write.csv(res_dlpfc_b, "corrplot/dlpfc_b.csv")
write.csv(res_dlpfc_m, "corrplot/dlpfc_m.csv")
write.csv(res_dlpfc_s, "corrplot/dlpfc_s.csv")

#Corrplot
logfolds <- as.matrix(read.csv(file = "data/log2folds.csv", row.names="X"))
#head(logfolds)
jpeg(file="corrplot.jpeg")
corrplot(cor(logfolds))
dev.off()

class(plotHeatmap)

# ---------------- PAPER ---------------- #

# the paper approach:

# LRT
#dds <- DESeq(dds, test="LRT", reduced=~1)
#res <- results(dds)


jpeg(file="s_nacc_plot.jpeg")
hist(res$pvalue, col= "green")
dev.off()

#res_ancg_b$log2FoldChange

#Rename Columns
#Log2Fold Plot on all of them 
corrplot(res_ancg_b$log2FoldChange)


## variance stabilizing
vsd <- vst(dds, blind=FALSE)
