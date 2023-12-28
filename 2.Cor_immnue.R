
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringi)
library(GSVA)
library(limma)

Expression_data <- read.csv("Expression_data.csv",row.names=1)
immune_gene <- read.table("GeneList.txt", fill =TRUE, header = TRUE, sep = "\t")
immune_gene_list <- list()
for(pathway_name in names(table(immune_gene$Category))){
    immune_gene_list[[pathway_name]] <- immune_gene[which(immune_gene$Category %in% pathway_name),1]
}
immport_data <- gsva(expr=as.matrix(Expression_data), gset.idx.list=immune_gene_list, 
        method="ssgsea", kcdf="Gaussian", ssgsea.norm=TRUE, parallel.sz=64)
save(immport_data, file = "immport_data.RData")
cor_asymmetry <- function(input1=NULL, input2=NULL, cor_method="spearman", merge_df=FALSE){

  cor_matrix <- matrix(NA, nrow=dim(input1)[2], ncol=dim(input2)[2])
  rownames(cor_matrix) <- colnames(input1)
  colnames(cor_matrix) <- colnames(input2)

  pvalue_matrix <- matrix(NA, nrow=dim(input1)[2], ncol=dim(input2)[2])
  rownames(pvalue_matrix) <- colnames(input1)
  colnames(pvalue_matrix) <- colnames(input2)
  for(i in 1:dim(input1)[2]){
    for(j in 1:dim(input2)[2]){
      cor_test_res <- cor.test(input1[,i], input2[,j], method=cor_method)
      cor_matrix[i, j] <- cor_test_res$estimate
      pvalue_matrix[i, j] <- cor_test_res$p.value
    }
  }
  res_list <- list(correlation=cor_matrix, pvalue=pvalue_matrix)
  if(merge_df){
    cor_pval_merge <- lapply(res_list, function(input){
      res <- data.frame(var1=rownames(input), input)
      res <- convArrType(inputArr=res, keyCols="var1", valueCols=colnames(res)[-1], wide2long=TRUE)
      return(res)
    })
    cor_res <- cor_pval_merge[[1]]
    colnames(cor_res)[3] <- "corr"
    pval_res <- cor_pval_merge[[2]]
    colnames(pval_res)[3] <- "pval"
    cor_pval_merge <- merge(cor_res, pval_res, by=c("var1", "variable"))
    colnames(cor_pval_merge) <- c("var1", "var2", "corr", "pval")
    res_list[["cor_pval_merge"]] <- cor_pval_merge
  }
  return(res_list)
}
cor_data <- cor_asymmetry(input1=t(immport_data), input2=Rhythm_score[,1,drop = FALSE], cor_method="spearman")


