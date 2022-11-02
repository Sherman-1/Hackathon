install.packages(c("BiocManager","curl","purrr","XML"))

BiocManager::install("DESeq2", force = TRUE)

library(DESeq2)
