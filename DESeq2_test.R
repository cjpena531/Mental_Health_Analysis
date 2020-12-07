# code pulled from vignette: http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start

# install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("apeglm")

# ---------------- IMPORTS ---------------- #

# import packages
library("pasilla")
library("DESeq2")

# ---------------- DATA ---------------- #

setwd("~/Documents/Github/DESeq2_Tutorial/")

# get data
cts <- as.matrix(read.csv(file = "final_gene_counts.csv", row.names="refseq"))
#read.csv(pasCts,sep="\t",row.names="gene_id")
coldata <- read.csv("SraRunTable.csv", row.names="Run")
#head(coldata)
coldata <- coldata[,c("clinical_diagnosis","age_at_death", "Brain_pH", "post.mortem_interval")]
coldata$clinical_diagnosis <- factor(coldata$clinical_diagnosis)
coldata$age_at_death <- factor(coldata$age_at_death)
coldata$Brain_pH <- factor(coldata$Brain_pH)
coldata$post.mortem_interval <- factor(coldata$post.mortem_interval)

# look at the data
head(cts,2)
#coldata

# not in same order! 
#rownames(coldata) <- sub("fb", "", rownames(coldata))

# the same samples
all(rownames(coldata) %in% colnames(cts))
# but not the same order!
all(rownames(coldata) == colnames(cts))

# sort to be in the same order
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

# ---------------- DESeqDataSet ---------------- #

# create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ clinical_diagnosis)
dds

# set up metadata
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

# ---------------- FILTERING ---------------- #

# pre-filter
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# factors in R!
dds$clinical_diagnosis <- factor(dds$clinical_diagnosis, levels = c("Control", "Major Depression", "Bipolar Disorder", "Schizophrenia"))

# ---------------- DEA ---------------- #

# Differential expression analysis
?DESeq
# carries out: estimation of size factors, estimation of dispersion: neg. binomial GLM
dds <- DESeq(dds)
res <- results(dds)
res

# Log fold change
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

resOrdered <- resLFC[order(resLFC$pvalue),]
sum(resOrdered$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

# ---------------- PAPER ---------------- #

# the paper approach:

# LRT
dds <- DESeq(dds, test="LRT", reduced=~1)
res <- results(dds)

## variance stabilizing
vsd <- vst(dds, blind=FALSE)
