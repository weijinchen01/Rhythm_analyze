
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringi)
library(GSVA)
library(limma)

Expression_data <- read.csv("Expression_data.csv",row.names=1)

Rhythm_gene_list <- list(Rhythm_score= c("RORA", "PER3", "PER1", "CRY2", "RORB", "PER2", "NR1D2", 
"RORC", "CRY1", "ARNTL", "NPAS2", "CLOCK", "NR1D1", "ARNTL2"))
Rhythm_score <- gsva(expr=as.matrix(Expression_data), gset.idx.list=Rhythm_gene_list, 
        method="ssgsea", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=64)
Rhythm_score <- as.data.frame(t(Rhythm_score))
save(Rhythm_score, file = "Rhythm_score.RData")

Rhythm_score$Group <- "Score_high"
Rhythm_score$Group[Rhythm_score[,1] < median(Rhythm_score[,1])] <- "Score_low"

pdata_tmp <- Rhythm_score[,"Group",drop = FALSE]
contrast_name=c("Score_high", "Score_low")
design_diff <- matrix(0, nrow=dim(pdata_tmp)[1], ncol=length(contrast_name))
  rownames(design_diff) <- rownames(pdata_tmp)
  colnames(design_diff) <- contrast_name
  for(i in 1:length(contrast_name)){
    design_diff[which(as.character(pdata_tmp[,1])==contrast_name[i]), i] <- 1
  }
fit <- lmFit(Expression_data, design_diff)
contrast.matrix <- makeContrasts(contrasts=paste(contrast_name, collapse="-"), levels=contrast_name)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit3 <- eBayes(fit2) 
diff_genes_Rhythm <- topTable(fit3, coef=1, adjust.method="none", number=dim(Expression_data)[1])

gmtfile <- "MsigDB/version7.5/h.all.v7.5.1.symbols.gmt"
TERM2GENE_set <- read.gmt(gmtfile)
colnames(TERM2GENE_set) <- c("ont", "gene")
set.seed(1235)
gene_ordered <- diff_genes_Rhythm[, "logFC"]
names(gene_ordered) <- rownames(diff_genes_Rhythm)
gene_ordered <- sort(gene_ordered, decreasing=TRUE)
res <- GSEA(gene_ordered, exponent=1, pAdjustMethod="none", TERM2GENE=TERM2GENE_set, seed=FALSE, by="fgsea")

