rm(list=ls())

library(ggplot2)
library(heatmaply)
library(corrplot)
library(reshape2)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(data.table)
library(Rtsne)
library(sva)
library(preprocessCore)
library(factoextra)
library(MCPcounter)
library(dplyr)
library(tidyverse)
library(drake)

setwd("C:/Users/aflorescu/Molecular Partners AG/DEV_TP_ExVivo - Ana/Rscripts/Rprojects/Analysis-bulk-cohorts")

#====================================
#Functions to transpose count matrix and join to annotation


transpose.TPM <- function(df.TPM){
  rownames(df.TPM) <- df.TPM$Gene
  df.TPM <- df.TPM[,-1]
  df.t <- data.frame(SAMPLE_ID=colnames(df.TPM),t(df.TPM))
  #df.merged<- merge(df.annotation, df.t, by="SAMPLE_ID")
  return(df.t)
  
}

merge.file.annotation <- function(df.t,df.annotation){
  #rownames(df.TPM) <- df.TPM$Gene
  #df.TPM <- df.TPM[,-1]
  #df.t <- data.frame(SAMPLE_ID=colnames(df.TPM),t(df.TPM))
  df.merged<- merge(df.annotation, df.t, by="SAMPLE_ID")
  return(df.merged)
  
}


#============================================================================================================
#main workflow

plan <- drake_plan(
  
  df.TPM = fread("C:/data/POG/POG570_TPM_expression_log2-add-1-combat-by-protocol-qnorm.txt", header=TRUE, sep="\t", data.table = FALSE),
  df.annotation=read.table("C:/data/JoinedFiles/TCGA-POG-MCP-GSVA-11.11.2020.txt", header = TRUE, sep = "\t"),
  
  df.t =transpose.TPM(df.TPM),
  df.joined = merge.file.annotation(df.t, df.annotation),
  
  pca.results = prcomp(df.t[,-1]),
  plot1 =fviz_eig(pca.results),
  
  coords.pca = data.frame(SAMPLE_ID=rownames(pca.results$x),pca.results$x[,1:10]),
  annotation.pca = merge(df.annotation, coords.pca, by="SAMPLE_ID"),
   

  report = target(
    command = {
      rmarkdown::render(knitr_in("analysis-report.Rmd"))
      file_out("analysis-report.html")
      
    }
  )
)


vis_drake_graph(plan)

make(plan)
