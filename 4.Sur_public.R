
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringi)
library(GSVA)
library(limma)
library(survminer)
load("Data_exprs_list_public.RData")
load("Data_pdata_list_public.RData")
interest_gene <- "RORA"

for(GEO_ID in names(Data_exprs_list_public)){
  exprs_interest <- Data_exprs_list_public[[GEO_ID]][rownames(Data_exprs_list_public[[GEO_ID]]) %in% interest_gene,,drop = FALSE]
  exprs_interest <- as.data.frame(t(exprs_interest))
  exprs_interest <- data.frame(barcode=rownames(exprs_interest), exprs_interest)
  exprs_clinical_merged <- merge(exprs_interest, Data_pdata_list_public[[GEO_ID]], by.x="barcode", by.y="GEO_ID")
  surv_res_list_gene <- list()
  for(gene in colnames(exprs_interest)[-1]){
    if("OS" %in% colnames(exprs_clinical_merged)){
        sur_col <- c("OS", "OS_time")
    }else{
        sur_col <- c("PFS", "PFS_time")
    }
    exprs_clinical_merged[,sur_col[1]] <- as.numeric(exprs_clinical_merged[,sur_col[1]])
    exprs_clinical_merged[,sur_col[2]] <- as.numeric(exprs_clinical_merged[,sur_col[2]])
    var_cols=gene
    labels=c("Low", "High")
    res.cut <- surv_cutpoint(data=exprs_clinical_merged, time=sur_col[2], event=sur_col[1], variables=gene, minprop=0.15, progressbar = TRUE)
    res.cat <- surv_categorize(res.cut, labels=labels)
    res.cat[, var_cols] <- factor(as.character(res.cat[, var_cols]), levels=labels)
    surv_res_list_gene[[gene]] <- SurvAnalysis_univar(inputArr=res.cat, out_prefix=NULL, 
    variable_col=gene, time_col=sur_col[2], status_col=sur_col[1], force_plot=TRUE, cutoff=0.05, 
    conf.int=FALSE, risk.table=TRUE, title=gene, legend_position="none", break.time.by=NULL, 
    width=7, height=9)
  }
  surv_res_list_GEO[[GEO_ID]] <- surv_res_list_gene
}


