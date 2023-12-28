
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringi)
library(GSVA)
library(limma)

Expression_data <- read.csv("Expression_data.csv",row.names=1)
clinical_data <- read.csv("Clinical_data.csv")

data_clinical_sub <- clinical_data[, c("Sample", "Response","PFS", "PFST.Month")]
RORA_HDAC3_DDX3X_group <- Expression_data[c("RORA", "HDAC3", "DDX3X", "CD274"),] %>% t() %>% 
                    as.data.frame() %>% 
                    rownames_to_column()
RORA_HDAC3_DDX3X_group_cli <- merge(RORA_HDAC3_DDX3X_group, data_clinical_sub, by.x="rowname", by.y = "Sample", all.x =TRUE)
RORA_HDAC3_DDX3X_group_cli <- RORA_HDAC3_DDX3X_group_cli %>% 
                    filter(Response != "<NA>")    
RORA_HDAC3_DDX3X_group_cli <- RORA_HDAC3_DDX3X_group_cli %>% mutate(RORA_HDAC3_DDX3X = DDX3X/HDAC3*RORA) %>% 
                    mutate(group = case_when(RORA_HDAC3_DDX3X > median(.$RORA_HDAC3_DDX3X) ~ "High", 
                            RORA_HDAC3_DDX3X <= median(.$RORA_HDAC3_DDX3X) ~ "Low",TRUE ~ NA))

library(plyr)
percent_matrix <- data.frame(table(RORA_HDAC3_DDX3X_group_cli$group,RORA_HDAC3_DDX3X_group_cli$Response))
percent_matrix<- ddply(percent_matrix,.(Var1),transform,percent=Freq/sum(Freq)*100) 
percent_matrix$label = paste0(sprintf("%.1f", percent_matrix$percent), "%")

predictor <- RORA_HDAC3_DDX3X_group_cli$group
response <- RORA_HDAC3_DDX3X_group_cli$Response
uni_predictor<-unique(RORA_HDAC3_DDX3X_group_cli$group) 
uni_response<-unique(RORA_HDAC3_DDX3X_group_cli$Response)
number_matrix<-matrix(0,length(uni_response),length(uni_predictor))
colnames(number_matrix)<-uni_predictor
row.names(number_matrix)<-uni_response
for(i in 1:length(uni_predictor))
{
  for(j in 1:length(uni_response))
  {
    number_matrix[j,i]<-length(intersect(which(response==uni_response[j]),which(predictor==uni_predictor[i])))
  }
}
chisq<-chisq.test(number_matrix,correct=F)
pvalue<-list(p_value=chisq$p.value)

var_cols="RORA_HDAC3_DDX3X"
labels=c("Low", "High")
res.cut <- surv_cutpoint(data=RORA_HDAC3_DDX3X_group_cli, time="PFST.Month", event="PFS", variables=var_cols, minprop=0.15, progressbar = TRUE)
res.cat <- surv_categorize(res.cut, labels=labels)
res.cat[, var_cols] <- factor(as.character(res.cat[, var_cols]), levels=labels)
SurvAnalysis_univar(inputArr=res.cat, out_prefix=NULL, 
    variable_col="RORA_HDAC3_DDX3X", time_col="PFST.Month", status_col="PFS", force_plot=TRUE, cutoff=0.05, 
    conf.int=FALSE, risk.table=TRUE, title=NULL, legend_position="none", break.time.by=NULL, 
    width=7, height=9)

cor_data <- cor_asymmetry(input1=RORA_HDAC3_DDX3X_group[,"RORA_HDAC3_DDX3X",drop = FALSE], 
input2=RORA_HDAC3_DDX3X_group[,"CD274",drop = FALSE], cor_method="pearson")
