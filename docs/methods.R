## ----setup, include=FALSE, echo=F, warning= F, message=F-----------------
knitr::opts_chunk$set(message = FALSE, 
                      warning = FALSE, 
                      error = FALSE, 
                      tidy = FALSE,
                      fig.align="center", 
                      dpi = 300, 
                      cache = FALSE,
                      progress=FALSE, 
                      quite = TRUE)

## ---- echo=FALSE, results="hide", warning=FALSE,message=FALSE------------
suppressPackageStartupMessages({
  if("TCGAWorkflow" %in% installed.packages()[,1]) {
    library(TCGAWorkflow)
  } else {
    devtools::load_all(".")
  }
})

## ---- eval=TRUE, message=FALSE, warning=FALSE, include=FALSE-------------
library(TCGAWorkflowData)
library(DT)
library(TCGAbiolinks)

## ---- eval=FALSE, include=TRUE-------------------------------------------
#  #----------------------------
#  # Obtaining DNA methylation
#  #----------------------------
#  # Samples
#  lgg.samples <- matched_met_exp("TCGA-LGG", n = 10)
#  gbm.samples <- matched_met_exp("TCGA-GBM", n = 10)
#  samples <- c(lgg.samples,gbm.samples)
#  
#  #-----------------------------------
#  # 1 - Methylation
#  # ----------------------------------
#  # For methylation it is quicker in this case to download the tar.gz file
#  # and get the samples we want instead of downloading files by files
#  query <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
#                    data.category = "DNA methylation",
#                    platform = "Illumina Human Methylation 450",
#                    legacy = TRUE,
#                    barcode = samples)
#  GDCdownload(query)
#  met <- GDCprepare(query, save = FALSE)
#  
#  # We will use only chr9 to make the example faster
#  met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
#  # This data is avaliable in the package
#  save("met.20.samples.GBM.LGG.chr9.rda")

## ----echo=TRUE, message=FALSE,warning=FALSE------------------------------
library(TCGAWorkflowData)
data("met.20.samples.GBM.LGG.chr9")
#----------------------------
# Mean methylation
#----------------------------
# Plot a barplot for the groups in the disease column in the
# summarizedExperiment object

# remove probes with NA (similar to na.omit)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))

TCGAvisualize_meanMethylation(met,
                              groupCol = "disease_type",
                              group.legend  = "Groups",
                              filename = "mean_lgg_gbm.png",
                              print.pvalue = TRUE)

