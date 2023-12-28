
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringi)
library(GSVA)
library(limma)

Expression_data <- read.csv("Expression_data.csv",row.names=1)
RORA_ex_data <- Expression_data["RORA", ] %>% t() %>% as.data.frame()
load("TIP.RData")
cor_data <- cor_asymmetry(input1=t(merge_cohort_TIP), input2=RORA_ex_data, cor_method="spearman")

