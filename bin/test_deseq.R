#### SMALL DESEQ TEST ####

args = commandArgs(trailingOnly=TRUE)

library(DESeq2)
library(MASS)


#### READ IN FILES AND CLEAN METADATA #####

path_to_counts <- args[1]

countsTable <- read.table(path_to_counts,header = TRUE)

rownames(countsTable) <- countsTable[,1]



# Stock des métadonnées 
metaData <- countsTable[,c(1:6)]

# Retrait des métadonnées de la table de comptage
countsTable <- countsTable[,c(-1,-2,-3,-4,-5,-6)]


## Retrait extention fichier du nom SRAI ##
colnames(countsTable) <- sub("\\.bam", "", colnames(countsTable))


# On arrange les colonnes par ordre de SRAI
countsTable = countsTable[, order(colnames(countsTable))]

# On vire les gènes nuls
non_null_genes = rowMeans(countsTable)>0
countsTable = countsTable[non_null_genes,]

#### CONDITIONS EXPERIMENTALES ####
conds<-factor(c("Mut","WT"))
colData=data.frame(condition=conds) # Format objet DESeq2




# Load model
dds<- DESeqDataSetFromMatrix(countsTable,colData,design = ~condition)
dds <- estimateSizeFactors(dds)

# Remove genes with less than 10 reads on normalised counts
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

res <- results(dds)


png("Figures/p_values.png", width = 800, height = 800)
ord_pval = res$pvalue[order(res$pvalue)]
hist(res$pvalue)
dev.off()
