# import packages
library("pasilla")
library("DESeq2")
library("apeglm")
library("corrplot")
library("dplyr")
library("RColorBrewer")
library("gplots")

deseq <- function(directory, gene_counts, col_data, control, disease){
  setwd(directory)
  cts <- as.matrix(read.csv(file=gene_counts, row.names="refseq"))
  coldata <- read.csv(col_data, row.names="Run")
  coldata <- coldata[,c("source_name","age_at_death", "Brain_pH", "post.mortem_interval")]
  coldata$source_name <- factor(coldata$source_name)
  coldata$age_at_death <- factor(coldata$age_at_death)
  coldata$Brain_pH <- factor(coldata$Brain_pH)
  coldata$post.mortem_interval <- factor(coldata$post.mortem_interval)
  cts <- cts[, rownames(coldata)]
  
  # create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~source_name)
  featureData <- data.frame(gene=rownames(cts))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds$source_name <- factor(dds$source_name, levels = c(control, disease))
  dds <- DESeq(dds)
  res <- results(dds)
  resname <- resultsNames(dds)
  resLFC <- lfcShrink(dds, coef=resname[2], type="apeglm")
  return(resLFC)

}