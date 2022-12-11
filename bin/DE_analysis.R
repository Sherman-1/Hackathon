#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(DESeq2)
library(MASS)


#### READ IN FILES AND CLEAN METADATA #####

path_to_counts <- args[1]

countsTable <- read.table(path_to_counts,header = TRUE)

rownames(countsTable) <- countsTable[,1]



# Stock and remove metadata 
metaData <- countsTable[,c(1:6)]
countsTable <- countsTable[,c(-1,-2,-3,-4,-5,-6)]


# Clean and order colnames
colnames(countsTable) <- sub("\\.bam", "", colnames(countsTable))
countsTable = countsTable[, order(colnames(countsTable))]

# Remove all null genes first
non_null_genes = rowMeans(countsTable)>0
countsTable = countsTable[non_null_genes,]

# Load experimental conditions 
conds<-factor(c("Mut","Mut","Mut","WT","WT","WT","WT","Mut"))
colData=data.frame(condition=conds) # Format objet DESeq2




#### DIFFERENTIAL EXPRESSION ANALYSIS ####

# Load model
dds <- DESeqDataSetFromMatrix(countsTable,colData,design = ~condition)
dds <- estimateSizeFactors(dds)

# Remove genes with less than 5 reasd in at least 2 o
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

# Run analysis
res <- results(dds)

#### PLOTS #### 

# Adjusted and standard p_values

png("Figures/p_values.png", width = 800, height = 800)
hist(res$pvalue)
dev.off()

png("Figures/adj_p_values.png", width = 800, height = 800)
hist(res$padj)
dev.off()



# Define volcano plot function

volcano <- function(alpha = 0.05, res) {
  
  colors <- densCols(res$log2FoldChange, -log10(res$padj))
  plot(res$log2FoldChange, -log10(res$padj), col=colors, panel.first=grid(),
       main="Volcano plot", xlab="log2(fold-change)", ylab="-log10 adjusted p-value",
       pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-1,1), col="brown")
  abline(h=-log10(alpha), col="brown")
  
  de_genes <- abs(res$log2FoldChange) > 1 & res$padj < alpha 
  text(res$log2FoldChange[de_genes],
       -log10(res$padj)[de_genes],
       lab=rownames(res)[de_genes ], cex=0.4)
  
  
}

# MAplot

png("Figures/MA_plot.png", width = 800, height = 800)
DESeq2::plotMA(res)
dev.off()

# Vplot
png("Figures/Volcano_plot.png", width = 800, height = 800)
vplot = volcano(alpha = 0.05, res)
dev.off()


# List all differentially expressed genes

genes = rownames(res)
de_genes <- abs(res$log2FoldChange) > 1 & res$padj < 0.05 
de_genes = genes[de_genes]
# Some genes do not fall into conditions and thus return NA
de_genes = de_genes[!is.na(de_genes)]

write.table(de_genes, "Figures/de_genes.txt")
