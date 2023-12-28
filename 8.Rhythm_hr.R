
library(ggplot2)
library(tidyverse)
library(plyr)
 library(survival)

Expression_data <- read.csv("Expression_data.csv",row.names=1)
clinical_data <- read.csv("Clinical_data.csv")
data_clinical_sub <- clinical_data[, c("Sample", "OS", "OST.Month")]
interest_gene <- c("RORA", "PER3", "PER1", "CRY2", "RORB", "PER2", "NR1D2", "RORC", "CRY1", "ARNTL", "NPAS2", "CLOCK", "NR1D1", "ARNTL2")
Rhythm_ex_data <- Expression_data[interest_gene, ] %>% t() %>% as.data.frame() %>% rownames_to_column()
Rhythm_data <- merge(Rhythm_ex_data, data_clinical_sub, by.x="rowname", by.y = "Sample", all.x =TRUE)
Rhythm_data[, "OS"] <- as.numeric(Rhythm_data[, "OS"])
Rhythm_data[, "OST.Month"] <- as.numeric(Rhythm_data[, "OST.Month"])
HR_res_list <- list()
for(gene in interest_gene){
    var_cols=gene
    labels=c("Low", "High")
    res.cut <- surv_cutpoint(data=Rhythm_data, time="OST.Month", event="OS", variables=gene, minprop=0.15, progressbar = TRUE)
    res.cat <- surv_categorize(res.cut, labels=labels)
    res.cat[, var_cols] <- factor(as.character(res.cat[, var_cols]), levels=labels)
    univ_formulas <- sapply(var_cols, function(x){
    as.formula(paste0("Surv(", "OST.Month", ", ", "OS", ")~", x))
  })
  univ_models <- lapply(univ_formulas, function(formula_used, inputArr=res.cat){
    coxph(formula=formula_used, data=inputArr)
  })
  univ_results <- lapply(univ_models, function(coxph_fit_res){
    HR <- round(exp(coxph_fit_res$coefficients), digits = 4)
    HR_CI <- round(exp(confint(coxph_fit_res)), digits = 4)
    pvalues <- summary(coxph_fit_res)$coefficients[,"Pr(>|z|)"]
    Desc_Str <- paste0(HR, " (", HR_CI[1], ", ", HR_CI[2], ")")
    res <- data.frame(Object=rownames(HR_CI), HR=HR, HR_CI, pvalue=pvalues, Description=Desc_Str)
    colnames(res) <- c("Object", "HR", "HR_CI_025", "HR_CI_975", "pvalue", "HR (95% CI)")
    return(res)
  })
  univ.result <- Reduce("rbind", univ_results)
  HR_res_list[[gene]] <- univ.result 
}
HR_res <- plyr::ldply(HR_res_list,data.frame)

