library("DESeq2")

cts <- as.matrix(read.csv(file = "data/filtered_gene_counts.csv", row.names="refseq"))
coldata <- read.csv("/datasets/srp073813/reference/SraRunTable.csv", row.names="Run")

coldata <- coldata[,c("clinical_diagnosis","age_at_death", "Brain_pH", "post.mortem_interval")]
coldata$clinical_diagnosis <- factor(coldata$clinical_diagnosis)
coldata$age_at_death <- factor(coldata$age_at_death)
coldata$Brain_pH <- factor(coldata$Brain_pH)
coldata$post.mortem_interval <- factor(coldata$post.mortem_interval)

cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                            design = ~ clinical_diagnosis)

# set up metadata
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$clinical_diagnosis <- factor(dds$clinical_diagnosis, levels = c("Control", "Major Depression", "Bipolar Disorder", "Schizophrenia"))

vsd <- vst(dds, blind=FALSE)

png(file="plots/pca.png",
width=600, height=600)
plotPCA(vsd, intgroup=c("clinical_diagnosis"))
dev.off()