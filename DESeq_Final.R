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

histogram <- function(var){
  name <- deparse(substitute(var))
  file_name <- paste(name, ".jpeg", sep="")
  jpeg(file="s_nacc_plot.jpeg")
  hist(var$pvalue, ylim = c(0, 4000))
  dev.off()
}

correlation_plot <- function(filename){
  logfolds <- as.matrix(read.csv(file = filename, row.names="X"))
  jpeg(file="corrplot.jpeg")
  corrplot(cor(logfolds))
  dev.off()
}

gene_matrix <- function(var) {
  var05 <- var[var$pvalue < 0.05,]
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  jpeg(file="heatmap.jpeg")
  heatmap.2((assay(var05)), col= hmcol, scale="row", Rowv = FALSE, Colv = TRUE)
  dev.off()
}

#DESeq2 Objects 
res_ancg_b <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_ancg_b.csv", "data/ancg_b.csv", "AnCg_Control", "AnCg_Bipolar Disorder")
res_ancg_m <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_ancg_m.csv", "data/ancg_m.csv", "AnCg_Control", "AnCg_Major Depression")
res_ancg_s <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_ancg_s.csv", "data/ancg_s.csv", "AnCg_Control", "AnCg_Schizophrenia")
res_nacc_b <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_nacc_b.csv", "data/nacc_b.csv", "nAcc_Control", "nAcc_Bipolar Disorder")
res_nacc_m <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_nacc_m.csv", "data/nacc_m.csv", "nAcc_Control", "nAcc_Major Depression")
res_nacc_s <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_nacc_s.csv", "data/nacc_s.csv", "nAcc_Control", "nAcc_Schizophrenia")
res_dlpfc_b <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_dlpfc_b.csv", "data/dlpfc_b.csv", "DLPFC_Control", "DLPFC_Bipolar Disorder")
res_dlpfc_m <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_dlpfc_m.csv", "data/dlpfc_m.csv", "DLPFC_Control", "DLPFC_Major Depression")
res_dlpfc_s <- deseq("~/Documents/Github/DESeq2_Tutorial/", "data/gc_dlpfc_s.csv", "data/dlpfc_s.csv", "DLPFC_Control", "DLPFC_Schizophrenia")
  
histogram(res_ancg_b)

#Make CSVs for Corrplot
write.csv(res_ancg_b, "corrplot/ancg_b.csv")
write.csv(res_ancg_m, "corrplot/ancg_m.csv")
write.csv(res_ancg_s, "corrplot/ancg_s.csv")
write.csv(res_nacc_b, "corrplot/nacc_b.csv")
write.csv(res_nacc_m, "corrplot/nacc_m.csv")
write.csv(res_nacc_s, "corrplot/nacc_s.csv")
write.csv(res_dlpfc_b, "corrplot/dlpfc_b.csv")
write.csv(res_dlpfc_m, "corrplot/dlpfc_m.csv")
write.csv(res_dlpfc_s, "corrplot/dlpfc_s.csv")

correlation_plot("data/log2folds.csv")

gene_matrix(res_ancg_b)

