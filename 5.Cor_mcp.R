
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringi)
library(GSVA)
library(limma)

Expression_data <- read.csv("Expression_data.csv",row.names=1)
RORA_ex_data <- Expression_data["RORA", ] %>% t() %>% as.data.frame()
library(IOBR)
mcp_res <- deconvo_tme(eset = Expression_data, method = "mcpcounter") %>% as.data.frame() %>% column_to_rownames("ID")
cor_data <- cor_asymmetry(input1=mcp_res, input2=RORA_ex_data, cor_method="spearman")



