## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ---- echo=FALSE, results="hide", warning=FALSE,message=FALSE------------
library(BiocStyle)
suppressPackageStartupMessages({
  if("TCGAWorkflow" %in% installed.packages()[,1]) {
    library(TCGAWorkflow)
  } else {
    devtools::load_all(".")
  }
})

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  deps <- c("pathview","clusterProfiler","ELMER", "DO.db","GO.db",
#            "ComplexHeatmap","EDASeq", "TCGAbiolinks","AnnotationHub",
#            "gaia","ChIPseeker","minet","BSgenome.Hsapiens.UCSC.hg19",
#            "MotifDb","MotIV", "rGADEM", "motifStack","RTCGAToolbox")
#  for(pkg in deps)  if (!pkg %in% installed.packages()) biocLite(pkg, dependencies = TRUE)
#  deps <- c("devtools","DT","pbapply","readr","circlize")
#  for(pkg in deps)  if (!pkg %in% installed.packages())  install.packages(pkg,dependencies = TRUE)
#  devtools::install_github("BioinformaticsFMRP/TCGAWorkflowData")
#  devtools::install_github("BioinformaticsFMRP/TCGAWorkflow", dependencies = TRUE)

## ---- eval=TRUE, message=FALSE, warning=FALSE, include=TRUE--------------
library(TCGAWorkflowData)
library(DT)

